from setuptools import setup

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='pyHETobs',
      version='0.2',
      description='Python Tool for calculating various HET parameters',
      long_description = readme(),
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Astronomy',
      ],
      keywords='HET Astronomy',
      url='https://github.com/indiajoe/pyHETobs',
      author='Joe Ninan',
      author_email='indiajoe@gmail.com',
      license='GPLv3+',
      packages=['pyHETobs'],
      install_requires=[
          'numpy',
          'astropy',
          'astroplan',
          'scipy',
          'matplotlib',
          'shapely'
      ],
      include_package_data=True,
      zip_safe=False)
