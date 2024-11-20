from setuptools import setup, find_packages

# 动态读取 requirements.txt 文件
def parse_requirements(filename):
    with open(filename, 'r') as file:
        return [line.strip() for line in file if line.strip() and not line.startswith('#')]

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
    install_requires=parse_requirements('./requirements.txt'),  # 从 requirements.txt 加载依赖
)
