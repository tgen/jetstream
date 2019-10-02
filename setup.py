import os
from setuptools import setup, find_packages


def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as fp:
        return fp.read()

setup(
    name='jetstream',
    version='1.5-b1',
    author='Ryan Richholt',
    author_email='ryan@tgen.org',
    url='https://github.com/tgen/jetstream',
    description='NGS analysis pipeline at TGen.',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    keywords='ngs pipeline automation',
    packages=find_packages(exclude=('test',)),
    include_package_data=True,
    python_requires='>=3',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Topic :: Utilities',
    ],
    install_requires=[
        'networkx',
        'pyyaml',
        'ulid-py',
        'jinja2',
        'filelock',
        'confuse'
    ],
    package_data={
        'jetstream': ['config_default.yaml']
    },
    entry_points={
        'console_scripts': [
            'jetstream=jetstream.cli:main',
        ],
    }
)
