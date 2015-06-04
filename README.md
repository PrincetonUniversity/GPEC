General Perturbed Equilibrium Codes
===================================

This is a repository for the suite of perturbed equilibrium package codes including IPEC, RPEC, GPEC, and PENTRC.

The primary purpose of these codes is to calcuated the nonaxisymmetric force balance in tokamak plasmas.


How to Contibute
-----------------

There is a `github repository <https://github.com/PrincetonUniversity/GPEC>`_ for this package to facilitate version control and collaborative contributions.

To become a contributor to the GPEC package, first go to the repository on `github <https://github.com/PrincetonUniversity/GPEC>`_ and fork using the button in the upper right. Copy the clone url from the box in the right side panel of your forked repository. Now create a local directory in your user area on the GA computers where you will make your changes. Navigate there and do::

   git clone <copied-url-from-github-fork> gpec
   cd gpec
   git checkout -b <new-branch>

and for convinience, add::

   git remote add github <copied-url-from-github-fork>

By convention, use your pppl username as the branch. Now you can freely edit files in the lib and bin directories. As you add files or make changes remember to::

   git add <new-file>
   git commit -a # Will open $EDITOR for message

When you have finished a significant addition/improvement push it using::

   git push --set-upstream github <new_branch>
   
When things are really working beautifully and you are ready to share your new developments with the community at large, go back on github, go to your forked repository, and submit a new pull request. Email nlogan@pppl.gov to have your code reviewed before being merged with the master. Continue to push to the same branch to update the pull request.

**Building the documentation**

The very first time, starting from the root of the repo::

    cd doc
    make init
    
With that set up, the following will build the documentation and push it to the proper place::

    make ghpages

Which is all you need do to update docs from there on out.

Documentation
-------------

Documentation of the GPEC package has been created using sphinx and is hosted by `PrincetonUniveristy GitHub Pages <http://princetonuniversity.github.io/GPEC/>`_.
