import os
from setuptools import setup, find_packages

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            rel_path = path.replace('jetstream/', '')
            paths.append(os.path.join(rel_path, filename))
    return paths

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


plugins = package_files('jetstream/plugins')
libs = package_files('jetstream/jslib')

setup(
    name = "jetstream",
    version = "0.1.0",
    author = "Ryan Richholt",
    author_email = "rrichholt@tgen.org",
    description = "NGS analysis pipeline at TGen.",
    long_description=read('README.md'),
    keywords = "ngs pipeline automation",
    packages=find_packages(),
    install_requires=['networkx', 'pydot'],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
    ],
    include_package_data=True,
    package_data={'': plugins},


    # Eventually this will be replaced with cli entry-points but for now it
    # helps me prototype quickly
    scripts=['examples/'+script for script in os.listdir('examples/')]

    # TODO develop cli for workflows
    # entry_points={
    #     'console_scripts': [
    #         'jetstream-run=jetstream.scripts.run:main',
    #         'jetstream-create=jetstream.scripts.run:main',
    #     ],
    # }

)
