#
# Copyright (C) 2020 RFI
#
# Author: James Parkhurst
#
# This code is distributed under the Apache license, a copy of
# which is included in the root directory of this package.
#

from setuptools import setup, find_packages


def main():
    """
    Setup the package

    """
    setup(
        package_dir={"": "src"},
        packages=find_packages(where="src"),
        install_requires=[
            "numpy",
            "python-parakeet",
        ],
        entry_points={"console_scripts": ["whippet=whippet.command_line:main"]},
        include_package_data=True,
    )


if __name__ == "__main__":
    main()
