fix charge_update command
========================

Syntax
""""""

.. parsed-literal::

   fix ID group-ID charge_udpate potential eta keyword values ...

* ID = name of fix
* group-ID = name of group fix is applied to
* potential = electric potential in Volts

  .. parsed-literal::

    *symm(etry) on/off*
        turn on/off charge neutrality constraint
    *couple <group1> potential1>*
        group1 = group of electrode atoms
        potential1 = electric potential in Volts
    *write_mat* <filename>
        write elestance matrix to file
    *write_inv* <filename>
        write inverted matrix to file
    *read_mat* <filename>
        read elestance matrix from file
    *read_inv* <filename>
        read inverted matrix from file


Examples
""""""""

.. code-block:: LAMMPS

   fix fxconp electrodes charge_update 0.0 1.805 symm on
   fix fxconp bot charge_update -1.0 1.805 couple top 1.0 write_inv inv.csv symm on

Description
"""""""""""

    The fix charge_update implements the constant potential method (CPM)
    (:ref:`Siepmann <Siepmann>`, :ref:`Reed <Reed3>`).
    Charges of groups specified as group-ID and with the `couple` keyword are
    adapted to meet their respective potential at every time step.
    An arbitrary number of electrodes can be set but the respective groups may
    not overlap.
    Electrode charges have a Gaussian charge distribution with reciprocal width
    eta.
    The energy minimization is achieved via matrix inversion :ref:`(Wang)
    <Wang5>`.
    The fix necessitates the use of a long range solver that can provide the 
    matrix of electrode-electrode interactions and a vector of 
    electrode-electrolyte interactions. 
    The Kspace styles *ewald/conp* and *pppm/conp* :ref:`(Ahrens-Iwers) <Ahrens>` are 
    created specifically for this task.

    For systems with non-periodic boundaries in one or two directions dipole
    corrections are available with the :doc:`kspace_modify <kspace_modify>`.
    For ewald/conp a two-dimensional Ewald summation can be used by setting
    "slab ew2d":

.. code-block:: LAMMPS

   kspace_modify slab <slab_factor>
   kspace_modify wire <wire_factor>
   kspace_modify slab ew2d

.. warning::

   Atom positions of electrode particles have to be fixed at all times.

----------

.. _Siepmann:

**(Siepmann)** Siepmann and Strik, J. Chem. Phys. 102, 511 (1995).

.. _Reed3:

**(Reed)** Reed *et al.*, J. Chem. Phys. 126, 084704 (2007).

.. _Wang5:
**(Wang)** Wang *et al.*, J. Chem. Phys. 141, 184102 (2014).

.. _Ahrens:

**(Ahrens-Iwers)** Ahrens-Iwers and Meißner, J. Chem. Phys. *in print*
