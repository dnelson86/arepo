/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/fof/fof_fuzz.c
 * \date        05/2018
 * \brief       Sorts particles which are not part of groups.
 * \details     contains functions:
 *                static void particle2in(data_in * in, int i, int firstnode)
 *                static void out2particle(data_out * out, int i, int mode)
 *                static void kernel_local(void)
 *                static void kernel_imported(void)
 *                void fof_fuzz_nearest_group(void)
 *                static int fof_fuzz_nearest_group_evaluate(int target, int mode, int threadid)
 *                void fof_assign_groups_to_fuzz()
 *                void fof_get_group_fuzz_offsets(void)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gsl/gsl_math.h>
#include <inttypes.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "../domain/domain.h"
#include "../main/allvars.h"
#include "../main/proto.h"
#include "../subfind/subfind.h"
#include "fof.h"

#if defined(FOF) & defined(FOF_FUZZ_SORT_BY_NEAREST_GROUP)

/* maximum search radius in units of Glob_r200_max */
#define FUZZ_MAX_GROUP_DISTANCE 10.0

static double Glob_r200_max;
static char *Todo;
static MyFloat *Hsml;
static double *NearestDist;

static void fof_get_group_fuzz_offsets(void);
static int fof_fuzz_nearest_group_evaluate(int target, int mode, int threadid);

static unsigned char *flag_node_contains_grp_particle;

/*! \brief Local data structure for collecting particle/cell data that is sent
 *         to other processors if needed. Type called data_in and static
 *         pointers DataIn and DataGet needed by generic_comm_helpers2.
 */
typedef struct
{
  MyDouble Pos[3];
  MyFloat Hsml;

  int Firstnode;
} data_in;

static data_in *DataIn, *DataGet;

/*! \brief Routine that fills the relevant particle/cell data into the input
 *         structure defined above. Needed by generic_comm_helpers2.
 *
 *  \param[out] in Data structure to fill.
 *  \param[in] i Index of particle in P and SphP arrays.
 *  \param[in] firstnode First note of communication.
 *
 *  \return void
 */
static void particle2in(data_in *in, int i, int firstnode)
{
#ifdef CELL_CENTER_GRAVITY
  if(P[i].Type == 0)
    {
      in->Pos[0] = SphP[i].Center[0];
      in->Pos[1] = SphP[i].Center[1];
      in->Pos[2] = SphP[i].Center[2];
    }
  else
#endif /* #ifdef CELL_CENTER_GRAVITY */
    {
      in->Pos[0] = P[i].Pos[0];
      in->Pos[1] = P[i].Pos[1];
      in->Pos[2] = P[i].Pos[2];
    }
  in->Hsml = Hsml[i];

  in->Firstnode = firstnode;
}

/*! \brief Local data structure that holds results acquired on remote
 *         processors. Type called data_out and static pointers DataResult and
 *         DataOut needed by generic_comm_helpers2.
 */
typedef struct
{
  double NearestDist;
  int GroupNr;
} data_out;

static data_out *DataResult, *DataOut;

/*! \brief Routine to store or combine result data. Needed by
 *         generic_comm_helpers2.
 *
 *  \param[in] out Data to be moved to appropriate variables in global
 *             particle and cell data arrays (P, SphP,...)
 *  \param[in] i Index of particle in P and SphP arrays
 *  \param[in] mode Mode of function: local particles or information that was
 *             communicated from other tasks and has to be added locally?
 *
 *  \return void
 */
static void out2particle(data_out *out, int i, int mode)
{
  if(mode == MODE_LOCAL_PARTICLES) /* initial store */
    {
      NearestDist[i] = out->NearestDist;
      PS[i].GroupNr  = out->GroupNr;
    }
  else /* merge */
    {
      if(out->NearestDist < NearestDist[i])
        {
          NearestDist[i] = out->NearestDist;
          PS[i].GroupNr  = out->GroupNr;
        }
    }
}

#include "../utils/generic_comm_helpers2.h"

/*! \brief Routine that defines what to do with local particles.
 *
 *  Calls the *_evaluate function in MODE_LOCAL_PARTICLES.
 *
 *  \return void
 */
static void kernel_local(void)
{
  int i;

  {
    int j, threadid = get_thread_num();

    for(j = 0; j < NTask; j++)
      Thread[threadid].Exportflag[j] = -1;

    while(1)
      {
        if(Thread[threadid].ExportSpace < MinSpace)
          break;

        i = NextParticle++;

        if(i >= NumPart)
          break;

        if(Todo[i])
          fof_fuzz_nearest_group_evaluate(i, MODE_LOCAL_PARTICLES, threadid);
      }
  }
}

/*! \brief Routine that defines what to do with imported particles.
 *
 *  Calls the *_evaluate function in MODE_IMPORTED_PARTICLES.
 *
 *  \return void
 */
static void kernel_imported(void)
{
  /* now do the particles that were sent to us */
  int i, cnt = 0;
  {
    int threadid = get_thread_num();

    while(1)
      {
        i = cnt++;

        if(i >= Nimport)
          break;

        fof_fuzz_nearest_group_evaluate(i, MODE_IMPORTED_PARTICLES, threadid);
      }
  }
}

/*! \brief Finds fuzzy particles to nearest group.
 *
 *  \return void
 */
void fof_fuzz_nearest_group(void)
{
  int i, npleft, iter = 0;
  long long ntot;

  CPU_Step[CPU_MISC] += measure_time();

  mpi_printf("FOF/FUZZ: Attaching fuzz...\n");

  Todo        = (char *)mymalloc("Todo", NumPart * sizeof(char));
  NearestDist = (double *)mymalloc("NearestDist", NumPart * sizeof(double));
  Hsml        = (MyFloat *)mymalloc("Hsml", NumPart * sizeof(MyFloat));

  for(i = 0; i < NumPart; i++)
    if(PS[i].GroupNr == 0)
      {
        NearestDist[i] = MAX_REAL_NUMBER;
        Hsml[i]        = All.ForceSoftening[P[i].SofteningType];
        Todo[i]        = 1;
      }
    else
      Todo[i] = 0;

  generic_set_MaxNexport();

  /* we will repeat the whole thing for those black holes where we didn't find enough neighbours */
  do
    {
      generic_comm_pattern(NumPart, kernel_local, kernel_imported);

      /* do final operations on results */
      for(i = 0, npleft = 0; i < NumPart; i++)
        {
          if(Todo[i])
            {
              if(PS[i].GroupNr < 0 || Hsml[i] > FUZZ_MAX_GROUP_DISTANCE * Glob_r200_max)
                Todo[i] = 0; /* done */
              if(PS[i].GroupNr == 0)
                {
                  Hsml[i] *= 2; /* increase search radius */
                  npleft++;
                }
            }
        }

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
        {
          iter++;

          if(iter > 0)
            mpi_printf("FOF/FUZZ: fuzz nearest group iteration %3d: need to repeat for %llu particles.\n", iter,
                       (unsigned long long)ntot);

          if(iter > MAXITER)
            {
              terminate("failed to converge in neighbour iteration\n");
            }
        }
    }
  while(ntot > 0);

  myfree(Hsml);
  myfree(NearestDist);
  myfree(Todo);

  mpi_printf("FOF/FUZZ: attaching fuzz done.\n");
  CPU_Step[CPU_FOF] += measure_time();
}

/*! \brief Evaluate function to find fuzzy particles to nearest group.
 *
 *  \param[in] target Index of particle/cell.
 *  \param[in] mode Flag if it operates on local or imported data.
 *  \param[in] threadid ID of thread.
 *
 *  \return Cost, i.e. number of nodes that had to be opened.
 */
static int fof_fuzz_nearest_group_evaluate(int target, int mode, int threadid)
{
  int no, k, numnodes, *firstnode;
  double dx, dy, dz, r, r2;
  MyDouble *pos;
  MyFloat h;
  double xtmp, ytmp, ztmp;
  double nearestDist = MAX_REAL_NUMBER;
  int groupNr        = 0;
  data_in local, *target_data;
  data_out out;

  if(mode == MODE_LOCAL_PARTICLES)
    {
      particle2in(&local, target, 0);
      target_data = &local;

      numnodes  = 1;
      firstnode = NULL;
    }
  else
    {
      target_data = &DataGet[target];

      generic_get_numnodes(target, &numnodes, &firstnode);
    }

  pos = target_data->Pos;
  h   = target_data->Hsml;

  for(k = 0; k < numnodes; k++)
    {
      if(mode == MODE_LOCAL_PARTICLES)
        {
          no = Tree_MaxPart; /* root node */
        }
      else
        {
          no = firstnode[k];
          no = Nodes[no].u.d.nextnode; /* open it */
        }

      while(no >= 0)
        {
          if(no < Tree_MaxPart) /* single particle */
            {
              dx = GRAVITY_NEAREST_X(Tree_Pos_list[3 * no + 0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Pos_list[3 * no + 1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Pos_list[3 * no + 2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;
              r  = sqrt(r2);
              if(r < nearestDist && r > 0)
                if(PS[no].GroupNr > 0)
                  {
                    nearestDist = r;
                    groupNr     = PS[no].GroupNr;
                  }
              no = Nextnode[no];
            }
          else if(no < Tree_MaxPart + Tree_MaxNodes) /* internal node */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                {
                  if(no <
                     Tree_FirstNonTopLevelNode) /* we reached a top-level node again, which means that we are done with the branch */
                    break;
                }

              struct NODE *current = &Nodes[no];
              int nocur            = no;
              no                   = current->u.d.sibling; /* in case the node can be discarded */

              if(nocur >= Tree_FirstNonTopLevelNode)
                if(flag_node_contains_grp_particle[nocur] == 0)
                  continue;

              double dist = h + 0.5 * current->len;
              dx          = FOF_NEAREST_LONG_X(current->center[0] - pos[0]);
              if(dx > dist)
                continue;
              dy = FOF_NEAREST_LONG_Y(current->center[1] - pos[1]);
              if(dy > dist)
                continue;
              dz = FOF_NEAREST_LONG_Z(current->center[2] - pos[2]);
              if(dz > dist)
                continue;

              /* now test against the minimal sphere enclosing everything */
              dist += FACT1 * current->len;
              if(dx * dx + dy * dy + dz * dz > dist * dist)
                continue;

              no = current->u.d.nextnode; /* ok, we need to open the node */
            }
          else if(no >= Tree_ImportedNodeOffset) /* point from imported nodelist */
            {
              int n = no - Tree_ImportedNodeOffset;

              dx = GRAVITY_NEAREST_X(Tree_Points[n].Pos[0] - pos[0]);
              dy = GRAVITY_NEAREST_Y(Tree_Points[n].Pos[1] - pos[1]);
              dz = GRAVITY_NEAREST_Z(Tree_Points[n].Pos[2] - pos[2]);

              r2 = dx * dx + dy * dy + dz * dz;
              r  = sqrt(r2);
              if(r < nearestDist && r > 0)
                if(PS[no].GroupNr > 0)
                  {
                    nearestDist = r;
                    groupNr     = PS[no].GroupNr;
                  }

              no = Nextnode[no - Tree_MaxNodes];
            }
          else /* pseudo particle */
            {
              if(mode == MODE_IMPORTED_PARTICLES)
                terminate("mode == MODE_IMPORTED_PARTICLES");

              if(target >= 0)
                tree_treefind_export_node_threads(no, target, threadid);

              no = Nextnode[no - Tree_MaxNodes];
            }
        }
    }

  out.GroupNr     = -groupNr;
  out.NearestDist = nearestDist;

  /* Now collect the result at the right place */
  if(mode == MODE_LOCAL_PARTICLES)
    out2particle(&out, target, MODE_LOCAL_PARTICLES);
  else
    DataResult[target] = out;

  return 0;
}

/*! \brief Main routine to attaching fuzzy particles to nearest group.
 *
 *  \return void
 */
void fof_assign_groups_to_fuzz()
{
  double t0, t1;
  int i;
  int flagged = 0;
  long long tot_flagged;
  double loc_r200_max = 0.0;

  mpi_printf("FOF/FUZZ: starting fuzz assignment...\n");

  /* get maximum R200 */
  for(i = 0; i < Ngroups; i++)
    if(Group[i].R_Crit200 > loc_r200_max)
      loc_r200_max = Group[i].R_Crit200;

  MPI_Allreduce(&loc_r200_max, &Glob_r200_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  mpi_printf("FOF/FUZZ: maximum R_Crit200=%g (maximum search radius=%g)\n", Glob_r200_max, FUZZ_MAX_GROUP_DISTANCE * Glob_r200_max);

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestOccupiedTimeBin);
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  /* create gravity tree */
  construct_forcetree(0, 0, 0, All.HighestOccupiedTimeBin); /* build tree for all particles */

#ifdef HIERARCHICAL_GRAVITY
  timebin_make_list_of_active_particles_up_to_timebin(&TimeBinsGravity, All.HighestActiveTimeBin);
#endif /* #ifdef HIERARCHICAL_GRAVITY */

  flag_node_contains_grp_particle =
      (unsigned char *)mymalloc("flag_node_contains_grp_particle", Tree_MaxNodes * sizeof(unsigned char));

  memset(flag_node_contains_grp_particle, 0, Tree_MaxNodes * sizeof(unsigned char));
  flag_node_contains_grp_particle -= Tree_MaxPart;

  /* tag node which contain group particle to speed up tree search */
  for(i = 0; i < NumPart; i++)
    {
      if(PS[i].GroupNr > 0)
        {
          int no = Father[i];
          while(no >= 0)
            {
              flag_node_contains_grp_particle[no] = 1;
              no                                  = Nodes[no].u.d.father;
            }
        }
    }

  for(i = Tree_MaxPart; i < Tree_MaxPart + Tree_MaxNodes; i++)
    if(flag_node_contains_grp_particle[i] == 1)
      flagged++;

  sumup_large_ints(1, &flagged, &tot_flagged);
  mpi_printf("FOF/FUZZ: total flagged group nodes = %lld\n", tot_flagged);

  /* do search */
  t0 = second();
  fof_fuzz_nearest_group();
  t1 = second();
  mpi_printf("FOF/FUZZ: finding fuzz group neighbors took = %g sec\n", timediff(t0, t1));

  /* free flags */
  myfree(flag_node_contains_grp_particle + Tree_MaxPart);

  /* free gravity tree */
  myfree(Father);
  myfree(Nextnode);
  myfree(Tree_Points);
  force_treefree();

  /* now set PS, assign proper values for sorting */
  for(i = 0; i < NumPart; i++)
    {
      if(PS[i].GroupNr < 0)
        PS[i].GroupNr = -PS[i].GroupNr;
      if(PS[i].GroupNr == 0)
        PS[i].GroupNr = TotNgroups + 1 + 1;
    }

  /* get the offsets */
  fof_get_group_fuzz_offsets();

  mpi_printf("FOF/FUZZ: fuzz assignment done.\n");
}

/*! \brief Get the offset of FoF group position into fuzz parts of snapshot.
 *
 *  \return void
 */
void fof_get_group_fuzz_offsets(void)
{
  int ngrp, n, j, type, i, off, ntype[NTYPES];
  int target, nexport, nimport;
  int *list_typeNumPart;
  long long *typeNumPartOffset;
  long long *keyBufexp, *keyBufimp, *list_min_keys, *list_max_keys;
  long long fuzz_min_key, fuzz_max_key, key;

  struct data_aux_sort *aux_sort = (struct data_aux_sort *)mymalloc("aux_sort", sizeof(struct data_aux_sort) * NumPart);

  /* global P[].Type=type offset */
  list_typeNumPart  = mymalloc("list_typeNumPart", NTask * sizeof(int));
  typeNumPartOffset = mymalloc("typeNumPartOffset", NTask * sizeof(long long));

  list_min_keys = mymalloc("list_min_keys", NTask * sizeof(long long));
  list_max_keys = mymalloc("list_max_keys", NTask * sizeof(long long));

  for(i = 0; i < NTYPES; i++)
    ntype[i] = 0;

  for(i = 0; i < NumPart; i++)
    {
      aux_sort[i].OriginTask  = ThisTask;
      aux_sort[i].OriginIndex = i;
      aux_sort[i].GrNr        = PS[i].GrNr;
#ifdef SUBFIND
      aux_sort[i].SubNr            = PS[i].SubNr;
      aux_sort[i].DM_BindingEnergy = PS[i].BindingEnergy;
#endif /* #ifdef SUBFIND */
      aux_sort[i].Type = P[i].Type;
#ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP
      aux_sort[i].key = PS[i].GroupNr;
#else  /* #ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP */
      aux_sort[i].ID = P[i].ID;
#endif /* #ifdef FOF_FUZZ_SORT_BY_NEAREST_GROUP #else */
      ntype[P[i].Type]++;
    }

  qsort(aux_sort, NumPart, sizeof(struct data_aux_sort), fof_compare_aux_sort_Type);

  for(i = 0, off = 0; i < NTYPES; off += ntype[i], i++)
    parallel_sort(aux_sort + off, ntype[i], sizeof(struct data_aux_sort), fof_compare_aux_sort_GrNr);

  for(type = 0; type < NTYPES; type++)
    {
      MPI_Allgather(&ntype[type], 1, MPI_INT, list_typeNumPart, 1, MPI_INT, MPI_COMM_WORLD);

      for(n = 1, typeNumPartOffset[0] = 0; n < NTask; n++)
        typeNumPartOffset[n] = typeNumPartOffset[n - 1] + list_typeNumPart[n - 1];

      /* min/max key of local fuzz (aux_sort fuzz is key-sorted) */
      for(i = 0, fuzz_min_key = 0, fuzz_max_key = 0; i < NumPart; i++)
        {
          if(aux_sort[i].Type != type)
            continue;

          /* fuzz */
          if(aux_sort[i].GrNr == TotNgroups + 1)
            {
              if(fuzz_min_key == 0)
                fuzz_min_key = aux_sort[i].key;
              fuzz_max_key = aux_sort[i].key;
            }
        }

      MPI_Allgather(&fuzz_min_key, sizeof(long long), MPI_BYTE, list_min_keys, sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);
      MPI_Allgather(&fuzz_max_key, sizeof(long long), MPI_BYTE, list_max_keys, sizeof(long long), MPI_BYTE, MPI_COMM_WORLD);

      for(n = 0; n < NTask; n++)
        Send_count[n] = 0;

      for(i = 0; i < Ngroups; i++)
        {
          key = Group[i].GrNr + 1;
          for(target = 0; target < NTask; target++)
            if(key >= list_min_keys[target] && key <= list_max_keys[target])
              {
                Send_count[target]++;
                break;
              }
        }

      MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

      for(j = 0, nimport = 0, nexport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
        {
          nexport += Send_count[j];
          nimport += Recv_count[j];

          if(j > 0)
            {
              Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
              Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
            }
        }

      keyBufexp = (long long *)mymalloc("keyBufexp", nexport * sizeof(long long));
      keyBufimp = (long long *)mymalloc("keyBufimp", nimport * sizeof(long long));

      for(n = 0; n < NTask; n++)
        Send_count[n] = 0;

      for(i = 0; i < Ngroups; i++)
        {
          key = Group[i].GrNr + 1;
          for(target = 0; target < NTask; target++)
            if(key >= list_min_keys[target] && key <= list_max_keys[target])
              {
                keyBufexp[Send_offset[target] + Send_count[target]] = key;
                Send_count[target]++;
                break;
              }
        }

      for(n = 0; n < NTask; n++)
        {
          for(j = Send_offset[n]; j < Send_offset[n] + Send_count[n]; j++)
            {
              if(keyBufexp[j] < list_min_keys[n] || keyBufexp[j] > list_max_keys[n])
                terminate("FOF_FUZZ_SORT_BY_NEAREST_GROUP: keyBufexp broken %llu %llu %llu", list_min_keys[n], keyBufexp[j],
                          list_max_keys[n]);
            }
        }

      /* exchange */
      for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          target = ThisTask ^ ngrp;

          if(target < NTask)
            {
              if(Send_count[target] > 0 || Recv_count[target] > 0)
                {
                  MPI_Sendrecv(&keyBufexp[Send_offset[target]], Send_count[target] * sizeof(long long), MPI_BYTE, target, TAG_KEY,
                               &keyBufimp[Recv_offset[target]], Recv_count[target] * sizeof(long long), MPI_BYTE, target, TAG_KEY,
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      /* closest fuzz for imported keys */
      for(j = 0, i = 0, n = 0; j < nimport; j++)
        {
          while(i < NumPart)
            {
              if(aux_sort[i].Type == type)
                {
                  if((keyBufimp[j] == aux_sort[i].key) && (aux_sort[i].GrNr == TotNgroups + 1))
                    {
                      keyBufimp[j] = n;
                      break;
                    }
                  if((keyBufimp[j] < aux_sort[i].key) && (aux_sort[i].GrNr == TotNgroups + 1))
                    {
                      keyBufimp[j] = 0;
                      break;
                    }
                  n++;
                }
              i++;
            }
        }

      /* add offset */
      for(j = 0; j < nimport; j++)
        if(keyBufimp[j] > 0)
          keyBufimp[j] += typeNumPartOffset[ThisTask];

      /* exchange */
      for(ngrp = 0; ngrp < (1 << PTask); ngrp++)
        {
          target = ThisTask ^ ngrp;

          if(target < NTask)
            {
              if(Send_count[target] > 0 || Recv_count[target] > 0)
                {
                  MPI_Sendrecv(&keyBufimp[Recv_offset[target]], Recv_count[target] * sizeof(long long), MPI_BYTE, target, TAG_KEY,
                               &keyBufexp[Send_offset[target]], Send_count[target] * sizeof(long long), MPI_BYTE, target, TAG_KEY,
                               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }
        }

      for(n = 0; n < NTask; n++)
        Send_count[n] = 0;

      for(i = 0; i < Ngroups; i++)
        {
          key = Group[i].GrNr + 1;
          for(target = 0; target < NTask; target++)
            if(key >= list_min_keys[target] && key <= list_max_keys[target])
              {
                Group[i].FuzzOffsetType[type] = keyBufexp[Send_offset[target] + Send_count[target]];
                Send_count[target]++;
                break;
              }
          if(target == NTask)
            Group[i].FuzzOffsetType[type] = 0; /* not found */
        }

      myfree(keyBufimp);
      myfree(keyBufexp);
    }
  myfree(list_max_keys);
  myfree(list_min_keys);
  myfree(typeNumPartOffset);
  myfree(list_typeNumPart);
  myfree(aux_sort);
}

#endif /* #if defined(FOF) & defined(FOF_FUZZ_SORT_BY_NEAREST_GROUP) */
