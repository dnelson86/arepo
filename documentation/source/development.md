Code development
************************


We strongly encourage further development of the code by other people. The idea
with this public version is to provide a well tested stable version of AREPO to 
the community as a basis of individual model development. Changes to the code can 
then be made in two ways: bugfixes and code extensions. The former is organized in 
the issue-tracking system of the repository, while for the 
second one, the developers should be contacted.

Issue reporrting
================

Problems with the code will in generally be reported to the issue tracker of the repository.
Therefore, if a problem occurs, the first thing to check is whether there already exists
an open issue, i.e. whether this is a known problem. If not, there are two ways to create 
a new issue.

 * The issue tracking system of the gitlab repository requires log-in to ``gitlab.mpcdf.mpg.de``.
   In general, AREPO users will not have this access, however, active developers may request
   a guest account there, and can then create own issues. These issues need to be specific
   and reproducible, or directly point to the problematic part of the code.
 * AREPO users without such an account should contact the authors in case of problems with the
   code. Also in this case examples that reproduce the issue greatly help and accellerate
   the problem-solving process.


Code extensions 
===============

We welcome code extensions by other people, however with certain requirements 
for the implemented modules. 

 * The modules are under the same terms of use as the AREPO code, i.e. a GPLv3 license. 
   All modules in this version of AREPO are free to use for everyone.
 * The number of interfaces to the main code should be as few as possible. 
 * The module should be of general interest. This implies in particular that 
   the implementation needs to work in a fully scalable fashion in massively-parallel 
   runs, with a minimum of shared/dublicated data.
 * The module should come with some additional examples for verification of its 
   functionality.

Developers interested to share their module with the community as a part of 
AREPO should contact the authors.


Major code updates
==================

Attached a list of important bugfixes and additions with the date and id of the
commit

+------------------------------------------------------------------+-----------------------+------------+
| **Description**                                                  | **date (dd.mm.yyyy)** | **commit** |
+==================================================================+=======================+============+
| Public version complete                                          | 20.03.2019            |            |
+------------------------------------------------------------------+-----------------------+------------+

