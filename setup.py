from setuptools import setup, find_packages

setup(
    name= 'MetaEVD68',
    version= '0.1',
    packages= find_packages(),
    include_packages_data= True,
    install_requires= [
        'click',
        'PyYaml',
        'pylatex',
        'FastQC',
        'Trim Galore!',
        'STAR',
        'Samtools',
        'Picard tools',
        'Bam-readcount',
        'Bcftools',
        'Bedtools',
        'Kraken'
    ],
    entry_points='''
        [console_scripts]
        MetaEVD68=MetaEVD68.main:cli
    ''',






)
