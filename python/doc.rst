Very short introduction to Python
=================================

Python is interpreted, dynamic-typed language.

..

   **Task 1.** Start Python interpreter by typing ``python`` to shell and
   type in *Hello world!*

   .. code-block:: python

      >>> # Let's try Hello world! This is a comment.
      >>> print "Hello world!"
      Hello world!

Basic datatypes
---------------

..

   **Task 2.** Try using Python as an interactive calculator. Play with basic
   arithmetic operations; check the distinction of ``int`` and ``float``
   datatypes. Complex numbers are written as ``1.0 + 2.0j``. Power is written
   using double asterisk. Elementary functions are available in ``math``
   module:

   .. code-block:: python

      >>> import math
      >>> math.sqrt(4.0)
      2.0

      >>> from math import cos
      >>> cos(0.0)
      1.0

      >>> # Import everything from math
      >>> # Generally not a preferred way - obfuscates the namespace
      >>> from math import *
      >>> sin(pi)
      1.2246467991473532e-16

      >>> # Let's learn how to fly
      >>> import antigravity

Container datatypes
-------------------

.. code-block:: python

      >>> # Lists are ordered collection of objects
      >>> philoshophers = ['Cimrman', 'Landau', 'Lifschitz', 'Souček']

      >>> # First item has index 0!
      >>> philosophers[0]
      'Cimrman'

      >>> # Items can accessed backwards by negative indices
      >>> print philosophers[-1]
      Souček

      >>> # Slice ending index is one item farther
      >>> philosophers[1:3]
      ['Landau', 'Lifschitz']

      >>> # Strides
      >>> philosophers[1:4:2]
      ['Landau', 'Sou\xc4\x8dek']

      >>> # Any index can be ommited
      >>> philosophers[:]
      ['Cimrman', 'Landau', 'Lifschitz', 'Sou\xc4\x8dek']
      >>> philosophers[::2]
      ['Cimrman', 'Lifschitz']

      >>> # Dictionary is indexed by any (hashable) object
      >>> glass_volume = {'wine': 0.2,
      ...                 'beer': 0.5,
      ...                 'slivovice': 0.05}
      >>> ethanol_concentation = {'wine': 0.1,
      ...                         'beer': 0.05,
      ...                         'slivovice': 0.52}
      >>> print "One glass of wine has", \
      ...    str(glass_volume['wine']*ethanol_concentration['wine']), \
      ...    "litres of ethanol."
      One glass of wine has 0.02 litres of ethanol.


Flow control and functions
--------------------------

.. code-block:: python

      >>> # Blocks are defined by indentation
      >>> # Use consistently spaces; don't mix with tabs - danger
      >>> for i in range(10, -1, -1):
      ...     print i, 'green', {True: 'bottle', False: 'bottles'}[i==1], \
      ...         'hanging on the wall'
      ...

      >>> for os in ['Windows', 'Linux', 'Apples MacOS X', 'BSD']:
      ...     if 's' in os:
      ...         print os, 'sucks'

      >>> a = [3, 7, 666, 42, 616]
      >>> # Find divisible by 3
      >>> a_3 = []
      >>> for n in a:
      ...     if n%3 == 0:
      ...         a_3.append(n)

      >>> # The same can be achieved by list comprehension
      >>> a_3 = [n for n in a if n%3 == 0]

      >>> def heaviside(x):
      ...     if x > 0.0:
      ...         y = 1.0
      ...     elif x < 0.0:
      ...         y = -1.0
      ...     else:
      ...         y = 0.0
      ...     return y

..

   **Task 3.** Exploiting ``glass_volume`` and ``ethanol_concentation``
   variables defined above write function taking dictionary with keys of
   beverage type and values of number of glasses drunk and returning total
   volume of alcohol drunk.

What is variable, mutabulity and imutability
---------------------------------------------

Every **variable** in Python **is just a name for an object**. (Remember,
evything in Python is object.) Understanding semantics of assignmenet operator
is crucial thing! Consider following snippet

.. code-block:: python

   >>> a = 42
   >>> b = a
   >>> b
   42
   >>> b = 666
   >>> a
   42

In this example ``a`` is a name for the integer object (with value 42). On the
second line name ``b`` was bound to the same object. Then name ``b`` was bound
to the other int object (with value 666). This cannot change the value of the
original object (which ``a`` bounds to). This holds for object of any type.

**The statement** ``name = object`` **causes that** ``name`` **afterwards has no
connection to the prior object it was referring to and the prior object is
not changed in any way.** (With the exception that original object may be
garbage-collected if referenced nowhere else.)

On the other hand there are of course ways to change (mutate) objects (which
are mutable). The basic numeric types like ``int``, ``float`` etc. are imutable.
Some container data types (for instance ``list``, ``dict``) and user-defined
objects (classes, see below) are mutable.

.. code-block:: python

   >>> drinks = ['beer', 'wine', 'wine', 'wine', 'cognac', 'wine']
   >>> drinks[0] = ['tea']
   >>> drinks
   ['tea', 'wine', 'wine', 'wine', 'cognac', 'wine']
   >>> drinks.append('last small congnac') # mutating object
   >>> drinks
   ['tea', 'wine', 'wine', 'wine', 'cognac', 'wine', 'last small congnac']

   >>> id(drinks)
   140457968276328
   >>> drinks = drinks + ['last small cognac']
   >>> # drinks now has expected value but is is a new object - compare id
   >>> drinks
   ['tea', 'wine', 'wine', 'wine', 'cognac', 'wine', 'last small congnac', 'last
   small congnac', 'last small congnac']
   >>> id(drinks)
   140457968276184

   >>> # operator += may mutate (mutable) object
   >>> drinks += ['last small cognac']
   >>> drinks, id(drinks)
   (['tea', 'wine', 'wine', 'wine', 'cognac', 'wine', 'last small congnac', 'last
   small congnac', 'last small congnac', 'last small congnac'], 140457968276184)

   >>> # On the other hand imutable object cannot be mutated in any way
   >>> a = 42
   >>> id(a)
   41239880
   >>> a += 1
   >>> id(a)
   41239856

Concluding, assignment operator everytime bounds a name on lhs to an object on
rhs so that nothing is mutated. Member methods (like ``list.append`` in the
example above) may mutate a mutable object. Operators like ``+=`` may mutate a
mutable object. It is depending on implementation. In fact, line
``drinks += ['last small cognac']`` is interpreted as
``drinks.__iadd__(['last small cognac'])`` which mutates the object. On the
other hand ``a = 42; a += 1`` is interpreted as
``a = 42; a = a.__add__(1)`` because ``int`` object has not ``__iadd__``
method (as it is imutable and cannot be incremented in-place) so that
``__add__`` method returning a new ``int`` object is called.

.. todo::

   Introduction to classes is needed to understand DOLFIN code.
