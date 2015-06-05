General Perturbed Equilibrium Codes
===================================

This package contains a suite of nonaxisymmetric stability and perturbed equilibrium codes including DCON, IPEC, RPEC, GPEC, and PENTRC.

The primary purpose of this package is to calculated the nonaxisymmetric force balance in tokamak plasmas.

Note that the DCON code was developed by Allen Glasser at Los Alamos National Laboratory, while IPEC and PENTRC were developed and are maintained by the GPEC team at the Princeton Plasma Physics Laboratory.

Installation
-------------

Public installations should be downloaded and kept up to date from the releases on Github. The releases remain private for now, so only those with access to the repository will be able to install the package.


Development
---------------

There is a `Github repository <https://github.com/PrincetonUniversity/GPEC>`_ for this package to facilitate version control and collaborative contributions.

To become a contributor to the GPEC package, contact nlogan@pppl.gov to be added to the development team. Next, go the the Github repository to familiarize yourself with the documentation, branch structure, latest commits, etc. While there, go ahead and copy the HTTPS clone URL from the box in the right side panel. Finally, create a directory on your local machine, clone the repository, and checkout the appropriate branch in order to start contributing.

On the PPPL portal computers, navigate to a directory in which you want to do your development (something like /p/gpec/users/<username>/gpec) and use following commands::

   git clone <copied-url-from-github> develop
   cd develop
   git checkout -b develop
   git pull origin develop

This will create a directory 'develop' and with the developmental branch of the repository in it. For more on branches, see the section below. Now you are free to add/edit files in this branch. Be sure to use:: 

   git add <new-file>

when creating new files, so that they come under version control. Be sure to make regular commits::

   git commit -a # Will open $EDITOR for message

as you go. You can't go wrong with more commits. Especially when developing as a team.

When you have finished a significant addition/improvement push it to Github using::

   git push -u origin develop

so that others working on the same branch can stay in sync with your developments. Similarly, you should be regularly pulling from the branch in order to stay in sync with what others are doing. To do so, use::

   git pull origin develop

There are many great git tutorials on the web, and great documentation for all the commands. Use them!

**GIT Workflow**

All developers need to use Vincent Driessen's `GitFlow <http://nvie.com/posts/a-successful-git-branching-model>`_ workflow when editing the GPEC package. PLEASE READ THE ENTIRE POST. It is short, very thorough, and good for both experienced and new git users.

The highlights are,
  - There are two permanent branches: master and develop
  - The master branch is only updated for at release ready stages
  - New features should be developed in short-lived (days) branches coming off of and merging back to the development branch.
  
Specific instructions are given in the link above as to exactly how to branch and merge these various branches. For example, the --no-ff option should be used when merging in order to keep branch histories. Just follow the examples and you wont go wrong!


Documentation
-------------

Documentation of the GPEC package has been created using sphinx and is hosted by `PrincetonUniveristy GitHub Pages <http://princetonuniversity.github.io/GPEC/>`_.

The documentation for the entire package is contained within the top level docs directory, and uses rst files with syntax that should be straight forward to mimic. If creating a new .rst file be sure to include it in the top level organization of the index.rst file.

**Building the documentation**

The very first time, starting from the root of the repo::

    cd docs
    make init
    
This creates a new gpec-docs directory in the same directory that contains the local repository (not in the repository itself) and checks out a devoted documentation branch. With that set up, the following will build the documentation and push it to the proper place::

    make ghpages

This is now all you need to do in order to update docs for that local repo from there on out.
