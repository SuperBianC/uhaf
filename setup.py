from setuptools import setup, find_packages

setup(
    name='uhaf',  
    version='0.0.1',
    author='Haiyang Bian',
    author_email='253273104@qq.com',  
    description='Unified Hierarchical Annotation Framework for Single-cell Data',  
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',  
    url='https://github.com/SuperBianC/uhaf', 
    packages=find_packages(),
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.5',
    install_requires=[
        'numpy==2.1.3',
        'pandas==2.2.3',
        ],
)
