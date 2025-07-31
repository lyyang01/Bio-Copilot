#!/usr/bin/env python
try:
    import win32api
except:
    # not available on linux
    pass
import os
import shutil
import sys
from os.path import join
from pathlib import Path

import yaml
from docopt import docopt
from pkg_resources import DistributionNotFound, get_distribution

from .defaultlog import log
from .project import Project, subprocess_run


def main():
    """
    Tool to configure python projects

    Usage::
        autogen [-s|-d|-M|-m|-p]

    options::
        -h, --help     Display this help message
        -v, --version  Show version and exit

        -s --setup     Generate setup.py only (default)
        -d, --docs     Generate sphinx docs
        -M,--major     Publish major release
        -m,--minor     Publish minor release
        -p,--patch     Publish patch release

    """
    # process args
    try:
        version = get_distribution("autogen").version
    except DistributionNotFound:
        version = "not installed"
    args = docopt(main.__doc__, version=version)
    log.info(args)

    p = Project()

    # docs
    if args["--docs"]:
        p.create_docs()
        return

    # release = update version; recreate setup.py; release to git; release to pypi
    if args["--major"] or args["--minor"] or args["--patch"]:
        if subprocess_run("git status --porcelain --untracked-files no"):
            log.error("git commit all changes before release")
            return
        if args["--major"]:
            p.update_version(0)
        elif args["--minor"]:
            p.update_version(1)
        elif args["--patch"]:
            p.update_version(2)
        p.create_setup()
        p.release()
    else:
        p.create_setup()


if __name__ == "__main__":
    main()
