# SMA-UMAT
In this folder, there are two UMATs available for you to simulate the constituive behaviours of SMAs using finite element software Abaqus.

1. Finite_Strain_SMA_UMAT.for, this file is the UMAT to simulate SMAs undergoing large deformations. See reference paper for the detail. 
Xu, Lei, Theocharis Baxevanis, and D. C. Lagoudas. "A three-dimensional constitutive model for the martensitic transformation in polycrystalline shape memory alloys under large deformation." Smart Materials and Structures 28.7 (2019): 074004. https://doi.org/10.1088%2F1361-665x%2Fab1acb

2. SMA-TRIP.for, this file is the UMAT that simulates the constitutive behaviour of SMAs considering the transformation-induced plasticity(TRIP) and the two-way shape memory effect(TWSME), see the following paper as reference.
Lei Xu, Alexandros Solomou, Theocharis Baxevanis, Dimitris Lagoudas, Finite strain constitutive modeling for shape memory alloys considering transformation-induced plasticity and two-way shape memory effect, International Journal of Solids and Structures,2020. https://doi.org/10.1016/j.ijsolstr.2020.03.009.

A simple tutorial with Abaqus input files are provided. The tutorial is for the UMAT without TRIP. However you can still follow this simple tutorial to get the SMA-TRIP response by 1) filling in the additional material properties for TRIP UMAT 2) Specifying the TRIP UMAT as the running subroutine in Abaqus. If you have any questions at this moment, feel free to send emails to sdf007xulei@gmail.com for help. You are highly welcomed by the authors to use these open-resource codes for your non-commericial usage, and please cite the aforementioned papers accordingly.

Best regards,

Lei Xu, Ph.D.

March 13, 2021

----------------------------------

New Product Development Engineer.

Schlumberger Technology Corp.

Houston, TX  77054

Email: sdf007xulei@gmail.com
