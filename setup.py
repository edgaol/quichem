from setuptools import find_packages, setup

setup(
   name='quichem',
   version='0.0.1',
   description='Quick chemicistry tools',
   author='Eduardo Gabriel Amaral de Oliveira',
   author_email='eduardo.g.amaral1997@gmail.com',
   packages=find_packages(),
   install_requires=['pint', 'sympy', 'periodictable']
)

