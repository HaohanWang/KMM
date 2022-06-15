from distutils.core import setup

setup(
      name='KMM',
      version='0.99',
      author = "Haohan Wang",
      author_email='haohanw@cs.cmu.edu',
      url = "https://github.com/HaohanWang/KMM",
      description = ' Gene Set Priorization Guided by Regulatory Networks with p-values through Kernel Mixed Model',
      packages=['model', 'util'],
      scripts=['kmm.py'],
    )
