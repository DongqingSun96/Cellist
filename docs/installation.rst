.. highlight:: shell

.. role:: bash(code)
   :language: bash

Installation
------------




System requirements
>>>>>>>>>>>>>>>>>>>

* Linux/Unix
* Python >= 3.10


We recommend to create an independent conda environment for STRIDE. If users do not have conda, please install Miniconda first:
::
   
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh


Install the stable version
>>>>>>>>>>>>>>>>>>>>>>>>>>

Step 1 Prepare conda environment for STRIDE.
::::::::::::::::::::::::::::::::::::::::::::
:: 

   conda create -n cellist python=3.10
   conda activate cellist

Step 2 Install STRIDE package from :bash:`pypi`.
::::::::::::::::::::::::::::::::::::::::::::::::
::

   pip install -U cellist


Install the developing version
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Step 1 Prepare conda environment for STRIDE.
::::::::::::::::::::::::::::::::::::::::::::
:: 

   conda create -n cellist python=3.10
   conda activate cellist

Step 2 Download STRIDE package from github.
:::::::::::::::::::::::::::::::::::::::::::
::

   git clone https://github.com/wanglabtongji/Cellist.git

Step 3 Install dependencies of STRIDE.
::::::::::::::::::::::::::::::::::::::
::

   cd Cellist
   pip install -r requirements.txt

Step 4 Install STRIDE.
::::::::::::::::::::::
::
  
   python setup.py install