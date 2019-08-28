Code development
************************


Scientific software development is a crucial aspect of computational (astro-)physics. 
With the progess made over the past decades, numerical methods as well as implementations
have evolved and increased in complexity. The scope of simulation codes has 
increased accordingly, such that the question how to organize development becomes important 
to scientific progress.

This version of Arepo is intended as a basis, providing an efficient code structure 
and state of the art numerical techniques for gravity and hydrodynamics. The code is 
completely documented and should allow computational astrophysicists to develop their 
own specialized modules on top of it. Practically, this is best done by hosting an
own repository for the development, which is occasonally updated from the main 
repository to include latest bugfixes. For this workflow to work properly, it is
important to keep the code-base relatively static. For this reason

**the base version of the code is not meant to be extended in functionality.**

This means that only bugfixes as well as updates in the documentation and examples will
be included.

Issue reporting
===============

A code of the scope of Arepo will ineviably have a number of issues and so far undetected bugs.
Detection of all issues is very complicated and time-consuming work. This means in practice that 
we rely on users that actually report all bugs and issues they come across to improve the quality  
of the code. We therefore encourage users to report all issues they have, including things
that remain unclear after reading the documentation. This way we hope to constantly improve the
qualtiy of the code in all its aspects.
 
To organize this, we created a code support forum (www.arepo-code.org/forums/forum/arepo-forum)
which is publically visible. To create posts on this forum, a registration with admin-approval
is required. We encourage users to sign up under www.arepo-code.org/register and make 
use of the forum whenever they have problems. 

Issues reported in the forum will then be confirmed by one of the authors and added to the 
repository issue tracker, which serves as a to-do list for bug fixes and improvements to 
the code. Therefore, if a problem occurs, the first thing to check is whether there already exists
an open issue, i.e. whether this is a known and confirmed problem.

**Please use the support forum instead of contacting the authors by email.**

**We encourage experienced users to provide help and answer some of the open questions on this forum.**


Code extensions 
===============

Extensions should be hosted in separate repositories or branches of the main repository.
We highly welcome such extensions and encourage developers to make them publically 
available (under GNU GPL v3). While this can be done independently of the authors, 
we would encourage developers to inform the authors once there is a publically 
available module in order to have a list of available extensions on the code homepage. 

Some guidelines for code extensions and new modules:

 * Modularity
 
  * Minimize the number of changes in existing source code files to a few function calls 
    and structured variables within existing elements.
  * The source code of an own modue should largely be in separate (i.e. new) files.
  
 * Documentation
 
  * All parameter and config-options should be clearly explained in this documentation.
    Feel free to add an addional page to this documentation explaining how your module
    works.
  * Document what each function does in the source code.
  * The Template-Config.sh file should have a short explanation for each flag.
  
 * Verification and examples
 
  * If possible, create one or more addional examples illustrating and testing the module.


