======================
Introduction to PyTFTB
======================

About PyTFTB
------------

The PyTFTB project began as a Python implementation of the
`TFTB toolbox <http://tftb.nongnu.org/>`_ developed by François Auger,
Olivier Lemoine, Paulo Gonçalvès and Patrick Flandrin. While the Python
implementation (henceforth referred to as PyTFTB) and the MATLAB implementation
(henceforth referred to as TFTB) are similar in the core algorithms and the
basic code organization, the very nature of the Python programming language
has motivated a very different approach in architecture of PyTFTB (differences
between the two packages have been discussed in detail in the next section).
Thus, someone who is familiar with TFTB should find the PyTFTB API comfortably
within grasp, and someone who is beginning with PyTFTB should find it a fully
self contained library to use.

Comparison of TFTB and PyTFTB
-----------------------------

TFTB is broadly a collection of MATLAB functions and demos that demonstrate the use of these functions.
A detailed reference of these functions can be found `here <http://tftb.nongnu.org/refguide.pdf>`_. The fact that
Python is a general purpose programming language affords the users and the developers a lot of freedom, especially with
regard to code reuse and interfacing. The important differences in implementation are as follows:

1. PyTFTB makes heavy use of Python's object oriented design. This allows for code reuse and interfacing. Algorithms
that are very closely related to each other can inherit from thew same base class and reuse each others methods.

2. In TFTB, visualization of time frequency distributions is handled by dedicated functions like `tfrview` and `tfrqview`
whereas in PyTFTB, they are tightly coupled to the specific representation being computed.

3. PyTFTB is heavily dependent on the SciPy stack - especially the NumPy and the SciPy libraries. Whichever piece of
code can be delegated to these libraries is delegated to them.

Quick Start
-----------

.. toctree::
   :maxdepth: 2

   quickstart/intro_examples_1.rst
   quickstart/intro_examples_2.rst

