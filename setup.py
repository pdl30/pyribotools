import os
from setuptools import setup, find_packages

setup(name='pyribotools',
      version='0.0.1',
      #packages=find_packages(),
      description='pyribotools is a Python module to analyze RIBO-seq NGS data',
      author='Patrick Lombard',
      author_email='ptk.lmb55@gmail.com',
      packages=['pyribotools'],
      package_data={"pyribotools":['data/*']},
      scripts=['scripts/pyribo_align.py', 'scripts/pyribo_utr.py', 'scripts/pyribo_profile.py'],
      install_requires=['pysam', 'pybedtools'],
      license='GPLv3',
      platforms='any',
      classifiers=[
         'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
         'Development Status :: 3 - Alpha',
         'Programming Language :: Python :: 2.7',
         'Environment :: Console',
      ],
      long_description="""

pyribotools is a Python module to analyze RIBO-seq NGS data

 Contact
=============

If you have any questions or comments about pyribotools, please feel free to contact me via
eMail: ptk.lmb55@gmail.com

""",
    )
