from distutils.core import setup

setup(
    name = 'alignment', 
    packages = setuptools.find_packages(),
    version = '0.1.0', 
    license = 'GPL', 
    description = 'Protein sequences alignment using the Smith-Waterman algorithm', 
    author = 'Yunzhe Jiang', 
    author_email = 'yunzhe.jiang@yale.edu',
    install_requires=['numpy', 'pandas'],
    classifiers=[
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8', 
    'Programming Language :: Python :: 3.9', 
    'Programming Language :: Python :: 3.10', 
    ],
)
