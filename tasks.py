"""
invoke file for automating build/config stuff.

Example:

    invoke build
    invoke --list

Can be executed everywhere inside the abinit-panel directory, including build directories.

Use: `pip install invoke --user` to install invoke package.
"""
import os
#import sys
import platform

from contextlib import contextmanager
try:
    from invoke import task
except ImportError:
    raise ImportError("Cannot import invoke package. Use `pip install invoke --user`")


HERE = os.path.dirname(__file__)


@contextmanager
def cd(path):
    """
    A Fabric-inspired cd context that temporarily changes directory for
    performing some tasks, and returns to the original working directory
    afterwards. E.g.,

        with cd("/my/path/"):
            do_something()

    Args:
        path: Path to cd to.
    """
    # Taken from monty.os
    cwd = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(cwd)


@task
def build(ctx):
    """
    Build the book.
    """
    with cd(HERE):
        ctx.run("jb build abipy_book", pty=True)

@task
def build_all(ctx):
    """
    Build the book.
    """
    with cd(HERE):
        ctx.run("jb build abipy_book --all", pty=True)


# When debugging your book build, the following options can be helpful:
#
# jupyter-book build -W -n --keep-going mybookname/
# This will check for missing references (-n), turning them into errors (-W), but will still attempt to run the full build (--keep-going), so that you can see all errors in one run.
#
# You can also use -v or -vvv to increase verbosity.

#@task
#def build_no_cache(ctx):
#    """
#    Build the docker image locally (Do not use cache when building the image)
#    """
#    _before_building()
#    with cd(HERE):
#        ctx.run("docker build --no-cache -t abinit-panel:latest .", pty=True)


#@task
#def submodules(ctx):
#    """Update submodules."""
#    with cd(HERE):
#        # https://stackoverflow.com/questions/1030169/easy-way-to-pull-latest-of-all-git-submodules
#        ctx.run("git submodule update --remote --init", pty=True)
#        #ctx.run("git submodule update --remote", pty=True)
#        ctx.run("git submodule update --recursive --remote", pty=True)

@task
def linkcheck(ctx):
    """
    If you’d like to make sure that the links outside of your book are valid,
    run the Sphinx link checker with Jupyter Book. This will check each of your external links
    and ensure that they resolve.
    """
    with cd(HERE):
        ctx.run("jupyter-book build abipy_book/ --builder linkcheck", pty=True)


@task
def du(ctx):
    """
    """
    with cd(HERE):
        ctx.run("du -msh -I _build abipy_book/", pty=True)

@task
def md2nb(ctx, md_filename):
    """
    """
    ctx.run("jupytext {md_filename} --to ipynb" ,pty=True)

