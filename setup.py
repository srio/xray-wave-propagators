
#
# Memorandum: 
#
# Install from sources: 
#     git clone https://github.com/srio/xray-wave-propagators
#     cd xray-wave-propagators
#     python -m pip install -e . --no-deps --no-binary :all:
#
#
# Use (in python):
# from sajidpropagators.prop import *
#
# Upload to pypi (when uploading, increment the version number):
#     python setup.py register (only once, not longer needed)
#     python setup.py sdist
#     python setup.py upload


from setuptools import setup

setup(name='sajidpropagators',
      version='0.1',
      description='sajid propagators',
      author='Sajid Ali',
      author_email='sajidsyed2021@u.northwestern.edu',
      url='https://github.com/s-sajid-ali/xray-wave-propagators',
      packages=["sajidpropagators.prop",
                ],
      install_requires=[
                        'matplotlib',
                        'numpy',
                        'numexpr',
                       ]
)

