# Basic git workflows with IBAMR

## Overview
Changes to IBAMR's source code are managed with the git version control
system. IBAMR's developers work on GitHub with the https://www.github.com/IBAMR/
organization.

GitHub, as the name implies, is a tool for working with git repositories: we use
GitHub to track things that need to be fixed or improved in IBAMR with the
[issue tracker](https://www.github.com/IBAMR/IBAMR/issues).

## Using git and IBAMR
IBAMR has about a dozen regular contributors, who, at any time, are working on
dozens of new features and fixes. We use git to track these changes, handle
merging of new features, and keep a log of all commited changes to the library.

Changes to IBAMR are peer-reviewed: instead of simply pushing new features as
they are finished, all changes undergo a period when the core IBAMR developers
review the new code and check it for possible bugs, style problems, and
consistency with the rest of the library. Git facilitates this workflow with
*branches*: changes to IBAMR are first commited to branches, which are then
considered (via GitHub's pull request mechanism) for merging with the master
branch of GitHub. This is a simplified variant of what is usually called the
[gitflow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow)
model: the primary difference being that our develop and master branches are the
same. Some branches for IBAMR are actively developed by multiple users; these
branches are located in the IBAMR repository on GitHub itself. It is also
possible to develop branches outside of the IBAMR organization: this is
discussed below in the *Forking, Cloning, and making a simple fix* example.

## Examples

### Forking, Cloning, and making a simple fix
If you are new to git and GitHub we recommend starting with the
[spoon-knife](https://github.com/octocat/Spoon-Knife) (all its missing is a
fork) example to see how forks and clones work.

The easiest way to get started with IBAMR is to create a GitHub account and
*fork* the official (that is, https://www.github.com/IBAMR/IBAMR) IBAMR
repository by clicking on the Fork button in the top right. This will create a
separate copy of the *IBAMR/IBAMR* repository in your own GitHub account. After
doing this, you can *clone* the repository to your local computer in the same
way as was done for *Spoon-Knife*.

Once you have a local copy of IBAMR you can start making changes. With a
brand-new repository one usually needs to set their authorship information:

> git config --local user.name "John Doe"
> git config --local user.email "john.doe@lab.gov"

It is usually helpful to tell the git repository on your computer about a few
other git repositories: these are called *remotes*. The most important remote is
your fork on GitHub (which is, itself, another git repository). The convention
is to call this remote the *origin*. Similarly, the official IBAMR repository on
GitHub is, by convention, called the *upstream* remote. If you cloned your GitHub
fork then that fork should be set as your origin:

```
> git remote show origin
* remote origin
  Fetch URL: https://github.com/john.doe/ibamr
  Push  URL: https://github.com/john.doe/ibamr
  HEAD branch: master
  [...]
```
if this isn't set up, then try running

> git remote add origin https://github.com/john.doe/ibamr

Use the same command to set up upstream:

> git remote add upstream https://github.com/ibamr/ibamr

That should be all the configuration your local repository needs: at this point
your authorship information is set up and the repository is correctly linked
with your copy of IBAMR on GitHub and the official IBAMR repository.

Suppose that you wanted to add another function to
`ibtk/src/lagrangian/LEInteractor.cpp`. You should start by creating a branch
with a relevant name. If we added a function that calculated the number of
interaction points, for example, we should name it something like

> git checkout -b add-leinteractor-count-ib-points

Once we are on the branch and have made changes we can stage changes by
executing

> git add ibtk/src/lagrangian/LEInteractor.cpp

and then commit them with a brief but explanatory message:

> git commit -m "Added a function for counting IB points."

You can inspect the commit by viewing the commit log:

> git log

which will show the commit message, commit hash, time and date, and authorship
information. Once this is done, one can put the commit on GitHub by running

> git push origin add-leinteractor-count-ib-points

and, navigating to https://github.com/ibamr/ibamr, one can then open a pull
request by hitting the button `Compare & Pull Request`.

### Cleaning up history on a branch
Sometimes its necessary to take a sequence of commits on a branch (here
`feature-branch`) and reorder, reword, or recombine them in some way. The best
way to do this is to do an *interactive rebase*:

> git checkout feature-branch
>
> git rebase -i master

This does two things:

1. It puts the commits in `feature-branch` on top of those in `master`: if there
   are merge conflicts, then git will stop and ask you to fix them.
2. The `-i` for (interactive) flag displays, in your editor, the commits and
   their hashes. Each line contains one commit: one can reorder commits by
   reordering the lines. Changing `pick` in the leftmost column will alter what
   is done with the commit:
   - `reword`: change the commit message
   - `fixup`: combine this commit with the one above it. The commit message is
     removed.
   - `squash`: combine this commit with the one above it. The commit message is
     kept.
   - `edit`: Keep the commit, but pause the rebase afterwards so that one can
     make additional changes.

After rebasing it is necessary to force push the newly edited commit history for
`feature-branch` to the remote branch. This can be achieved using:

> git push -f origin feature-branch

The commit history for the remote branch should now reflect the edited local 
changes made to `feature-branch`.   

### Keeping a branch up to date
One can also use rebasing to keep a branch up to date with another branch:
suppose you have made lots of commits on `feature-branch` and other developers
have made commits on `master`. To get commits on master into `feature-branch`
you should run

> git checkout feature-branch
>
> git rebase master

This will 'replay' the commits in `feature-branch` on top of (that is, after all
of it's commits) `master`. Rebasing is preferable to merging since it keeps all
related commits together. Merging interleaves the commits on `master` and
`feature-branch`: this is undesirable since it makes it difficult to keep track
of which commits are unique to `feature-branch`.

If there are more commits in another repository that you want on the current
branch then you can run, for example,

> git pull origin feature-branch

to get commits from `origin`'s branch named `feature-branch` in the current
branch.

## Tips and Tricks
- Most installations of git include a simple branch and commit visualization
  tool named `gitk`; check it out!
- Using git from the command line is difficult. We recommend using a graphical
  interface like Git Cola, GitKraken, Tower, or (for emacs users) magit.

## Additional Material
There are lots of tutorials on the Internet that describe possible workflows
with git. We suggest poking around Stack Overflow for information on more
complex git commands. Another good resource are the video lectures put out by
Wolfgang Bangerth, which are available on YouTube in
[lecture 32.75](https://www.youtube.com/watch?v=kqb3aIakftA) and
[lecture 32.8](https://www.youtube.com/watch?v=kAqp2hhv-DU). The first video
gives an overview into how git works and the second discusses how git and GitHub
work together.
