from setuptools import setup

setup(name='correlate',
      version='1.0',
      description='Correlated Mutation Tool',
      url='https://github.com/andreubp/good_correlated.git',
      author='Andreu Bofill & Marina Reixachs',
      author_email='andreu.bofill@gmail.com',
      license='CMP',
      install_requires=[
        'numpy',
        'Biopython',
        'matplotlib',
        'plotly',
      ],
      scripts=['bin/correlate'],
      packages=['modules'],
      zip_safe=False)
