from setuptools import setup, find_packages

setup(
    name='svphaser',
    version='1.0.0',
    packages=find_packages(),
    install_requires=[
        'pysam',
        'pandas',
    ],
    entry_points={
        'console_scripts': [
            'svphaser=svphaser.cli:main'
        ]
    },
    author='Pranjul Mishra, Sachin Gadakh',
    author_email='pranjul.mishra@proton.me',
    description='SvPhaser: Structural Variant Phasing Tool for Phased BAM + Unphased SV VCF',
    url='https://github.com/SFGLab/Team5_Phased_SV_Analysis/SvPhaser',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
