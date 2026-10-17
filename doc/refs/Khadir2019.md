Optical occiety Journal of t A 

OPTICS, IMAGE SCIENCE, AND VISION 

# Quantitative model of the image of a radiating dipole through a microscope

SAMIRA KHADIR <sup>1</sup>,* PATRICK C CHAUMET <sup>1</sup> GUILLAUME BAFFOU <sup>1</sup> AND ANNE SENTENAC<sup>1</sup>,<sup>2</sup> 

<sup>1</sup>Institut Fresnel, Aix Marseille Univ, CNRS, Centrale Marseille, Marseille, France 

<sup>2</sup>e-mail: anne.sentenac@fresnel.fr 

*Corresponding author: samira.khadir@fresnel.fr 

Received 23 November 2018; revised 1 February 2019; accepted 1 February 2019; posted 4 February 2019 (Doc. ID 352169); published 8 March 2019 

In this paper, we introduce a formalism to determine the relationship between the full vectorial electric field existing at the object plane of a microscope and that existing at the image plane. The model is then used to quantitatively simulate, in both phase and intensity, the image of a radiating electric dipole placed either in a homogeneous medium or in the vicinity of a substrate. These simulations are compared with experimental measurements on single gold nanoparticles carried out by quadriwave lateral shearing interferometry. © 2019 Optica Society of America 

https://doi.org/10.1364/JOSAA.36.000478 

## 1. INTRODUCTION

The point spread function (PSF) of a microscope, defined as the image ofa point source, such as a nanoparticle, a fluorescent bead, or a molecule [1,2], is a key concept of image formation theory. Several approaches have been developed to calculate the theoretical PSF ofa microscope. The first models were based on a scalar diffraction theory, like in the work of Richards and Wolf in 1959 [3] and others [4,5]. These models were limited to unpolarized point sources imaged with low numerical aperture objectives. Thus, they did not provide complete information about the imaged nano-object [6]. To take into account the polarization state of the dipole and generalize the development to high numerical apertures, more accurate models based on vectorial diffraction theory have been proposed [7–9]. However, most of the above-mentioned formalisms do not use Cartesian coordinates, but polar coordinates, as a legacy of the time where the fast Fourier transform did not exist and where calculating a simple integral was simpler than performing twodimensional Fourier transforms. In addition, the calculations are rarely carried out quantitatively. They give the intensity at the image plane up to a given constant factor [9], which prevents its use in metrology applications. Yet, the latter is necessary for interpreting the results of quantitative phase imaging techniques, such as digital holographic microscopy, shearing interferometries, or polarization microscopy [10]. 

In this work, we introduce a full-vectorial theory to quantitatively calculate the image of a radiating electric dipole, in intensity and phase, through a standard microscope in a 4-f mounting. Our approach, using Cartesian coordinates and electromagnetic treatment, is close to that developed in [11] for Cassegrain objectives and microspectroscopy applications. An introductory section of this paper derives the transforma tion of a simple plane wave with a given incidence angl through a microscope. As any beam propagating in free space can be expanded into a sum ofplane waves, this first part allow us to draw the link between any vectorial field at the object focal plane of a microscope and that at the image focal plane. This approach is applied in a second section to determine the field radiated by an electric dipole in a homogeneous medium and in the presence of a substrate. Finally, results of simulations using this model are compared with experimental images o gold nanoparticles obtained by quadriwave lateral shearing interferometry (QLSI). The Matlab code corresponding to the model introduced in this work is provided (see Code 1, Ref. [12]). 

## 2. DESCRIPTION OF THE THEORETICAL MODEL

## A. Image of a Plane Wave through a Microscope

We consider a microscope with an objective of focal length $f _ { o }$ and numerical aperture NA, and a tube lens with focal length $\boldsymbol { f } _ { t }$ in a $4 f$ configuration satisfying the Abbe sine condition [13]. We call Oz the optical axis of the system, oriented from the objective to the tube lens (see Fig. 1), and M the magni fication factor defined by $M = - f _ { t } / \bar { f } _ { o } .$ M is a negative quantity as the image is inverted with respect to the object. The object focal plane is located at $z = z _ { o }$ while the image focal plane is placed at z z<sub>i</sub>. We call object space the homogeneous medium before the objective and we assume that it has a refractive index n. We call image space the homogeneous medium that lies between the tube lens and the sensor and assume that it has a refractive index of 1 as it usually consists of air. 

We first consider the action of the microscope on a monochromatic plane wave characterized by its wave vector ${ \bf k } =$ $( \mathbf { k } _ { | | } , \gamma )$ with ${ \bf k } _ { \parallel } = ( k _ { x } , k _ { y } ) , \gamma = \sqrt { n ^ { 2 } k _ { 0 } ^ { 2 } - k _ { \parallel } ^ { 2 } }$ , and $k _ { 0 } = 2 \pi / \lambda$ is the wavenumber in vacuum. In the object space, the incoming electric field reads 

$$
\mathbf {E} _ {\mathrm{ob}} (\mathbf {r} _ {\parallel}, z) = \mathbf {e} (\mathbf {k} _ {\parallel}) \exp [ i \mathbf {k} \cdot (\mathbf {r} - \mathbf {r} _ {o}) ],\tag{1}
$$

where $\mathbf { r } _ { o } = ( 0 , 0 , z _ { o } )$ is the object focal point, which sets the phase origin in the object space. The complex vectorial amplitude ${ \bf e } ( { \bf k } _ { \parallel } )$ of the plane wave can be decomposed on the polarization basis s, p defined as $\begin{array} { r } { \mathbf { s } = \frac { \mathbf { z } \times \mathbf { k } } { | \mathbf { z } \times \mathbf { k } | } } \end{array}$ <sup>and</sup> <sup>p</sup>  <sup>k×s</sup>k×s <sup>(see</sup> <sup>Fig.</sup> <sup>1).</sup> If $k _ { \parallel } \leq \mathrm { N A } k _ { 0 }$ , the microscope transforms the plane wave with wavevector k into a plane wave with wavevector $\mathbf { k } ^ { \prime }$ with $\mathbf { k } ^ { \prime } = ( \mathbf { k } _ { | | } ^ { \prime } , \boldsymbol { \gamma } ^ { \prime } )$ , where $\mathbf { k } _ { \| } ^ { \prime } = ( k _ { x } / M , k _ { y } / M )$ and $\gamma ^ { \prime } =$ $\sqrt { k _ { 0 } ^ { 2 } - k _ { \parallel } ^ { \prime 2 } }$ (the image space being air). This transformation preserves the s component of the field throughout its propagation and rotates the $\boldsymbol { \mathscr { P } }$ component (see Fig. 1). Thus, in the image space, i.e., after the tube lens, the field reads 

$$
\mathbf {E} _ {\mathrm{im}} (\mathbf {r} _ {\parallel}, z) = \frac {1}{M} \sqrt {\frac {\gamma}{\gamma^ {\prime}}} b (\mathbf {k} _ {\parallel}) \mathbf {e} ^ {\prime} (\mathbf {k} _ {\parallel}) \exp [ i \mathbf {k} ^ {\prime} \cdot (\mathbf {r} - \mathbf {r} _ {i}) ],\tag{2}
$$

where $\mathbf { r } _ { i } = ( 0 , 0 , z _ { i } )$ is the image focal point (the conjugate point of $\mathbf { r } _ { o } )$ , which sets the phase origin in the image space, 

$$
\mathbf {e} ^ {\prime} (\mathbf {k} _ {\parallel}) = [ (\mathbf {e} (\mathbf {k} _ {\parallel}) \cdot \mathbf {s}) \mathbf {s} + (\mathbf {e} (\mathbf {k} _ {\parallel}) \cdot \mathbf {p}) \mathbf {p} ^ {\prime} ],\tag{3}
$$

where $\begin{array} { r } { \mathbf { p } ^ { \prime } = \frac { \mathbf { k } ^ { \prime } \times \pmb { s } } { | \mathbf { k } ^ { \prime } \times \pmb { s } | } } \end{array}$ , and $b ( \mathbf { k } _ { \parallel } )$ is the cutofffunction ofthe microscope objective. In the case of a perfect microscope, this function is defined as $b ( \mathbf { k } _ { | | } ) = 1$ for $\begin{array} { r } { | \hat { \mathbf { k } } _ { \parallel } | < k _ { 0 } \mathrm { N A } , } \end{array}$ and 0 elsewhere. However, it should be noted that real microscope objectives exhibit an appodization effect that should be taken into account in the cutoff function for accurate predictions [14]. This effect depends on the microscope objective characteristics [15,16]. The leading factor $M ^ { - 1 } \sqrt { \frac { \gamma } { \gamma ^ { \prime } } }$ in Eq. (2) ensures energy conservation between the object space and the image space. Its expression is derived in Appendix A. 

Note that Eq. (3) is equivalent to the rotation of ${ \bf e } ( { \bf k } _ { | | } )$ : 

$$
\mathbf {e} ^ {\prime} (\mathbf {k} _ {\parallel}) = \mathbf {R} (\mathbf {k} _ {\parallel}) \mathbf {e} (\mathbf {k} _ {\parallel}),\tag{4}
$$

where R, the rotation matrix, is given in Appendix B. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-10-16/79454edb-07dd-45d2-b477-d4dc46f9cffc/2242edf0a1ce9eaf5ea9ea5281eaf8eeba85bc0d607f7f90456a795b9119c5d4.jpg)



Fig 1 Illustration of the effect of an optical microscope on a plane wave with a wavevector k.


## B. Link between the Image and Object Fields in a Microscope

Equations (1) and (2) are sufficient to thoroughly describe th image ofany object through the microscope, as any electromag netic beam propagating toward the positive z in the object space can be written as a sum of plane waves: 

$$
\mathbf {E} _ {\mathrm{ob}} (\mathbf {r}) = \iint \mathbf {e} (\mathbf {k} _ {\parallel}) \exp [ i \mathbf {k} \cdot (\mathbf {r} - \mathbf {r} _ {o}) ] \mathrm{d} \mathbf {k} _ {\parallel}.\tag{5}
$$

In this case, the field at any point in the image space (after the tube lens) reads 

$$
\mathbf {E} _ {\mathrm{im}} (\mathbf {r}) = \frac {1}{M} \iint \sqrt {\frac {\gamma}{\gamma^ {\prime}}} b (\mathbf {k} _ {\parallel}) \mathbf {e} ^ {\prime} (\mathbf {k} _ {\parallel}) \exp [ i \mathbf {k} ^ {\prime} \cdot (\mathbf {r} - \mathbf {r} _ {i}) ] \mathrm{d} \mathbf {k} _ {\parallel}.\tag{6}
$$

We now focus on the expression of the field at the image plane $[ \mathbf { r } - \mathbf { r } _ { i } = ( \mathbf { r } _ { \parallel } , 0 )$ in Eq. (6)] as it normally corresponds to the detector place. Since $\mathbf { k } _ { | | } = M \mathbf { k } _ { | | } ^ { \prime }$ and $\begin{array} { r } { \dot { { \bf d } } { \bf k } _ { \parallel } = { \bf \bar { \cal M } } ^ { 2 } \mathrm { d } { \bf k } _ { \parallel } ^ { \prime } , } \end{array}$ the electric field at the image plane reads 

$$
\mathbf {E} _ {\mathrm{im}} \left(\mathbf {r} _ {\parallel}, z _ {i}\right) = M \iint \sqrt {\frac {\gamma}{\gamma^ {\prime}}} b \left(M \mathbf {k} _ {\parallel} ^ {\prime}\right) \mathbf {e} ^ {\prime} \left(M \mathbf {k} _ {\parallel} ^ {\prime}\right) \exp \left(i \mathbf {k} _ {\parallel} ^ {\prime} \cdot \mathbf {r} _ {\parallel}\right) \mathrm{d} \mathbf {k} _ {\parallel} ^ {\prime}.\tag{7}
$$

Thus, the electric field at the image plane is the 2D Fourier transform of $\sqrt { \frac { \gamma } { \gamma ^ { \prime } } } b ( M \mathbf { k } _ { \parallel } ^ { \prime } ) \mathbf { e } ^ { \prime } ( M \mathbf { k } _ { \parallel } ^ { \prime } )$ , which is related to the elec tric field in the Fourier space of the microscope (i.e., the back focal plane of the objective). 

## C. Point Source in a Homogeneous Medium

We now consider the case of the image of a radiating dipole. The field radiated by a point source Q in a homogeneous space is the well-known solution of the equation [8] 

$$
\nabla \times (\nabla \times \mathbf {E} _ {Q}) - n ^ {2} k _ {0} ^ {2} \mathbf {E} _ {Q} = \mathbf {Q} \delta (\mathbf {r} - \mathbf {r} _ {Q}),\tag{8}
$$

which satisfies outgoing wave boundary conditions. The sourc Q relates to the electric dipole moment μ via 

$$
\mathbf {Q} = \frac {k _ {0} ^ {2}}{\varepsilon_ {0}} \boldsymbol {\mu}.\tag{9}
$$

$\mathbf { E } _ { Q }$ can be decomposed as a sum of plane waves using the $\tt W e y l$ expansion (i.e., the angular spectrum representation) [8,17]: 

$$
\mathbf {E} _ {Q} (\mathbf {r}) = \frac {i}{8 \pi^ {2}} \iint \frac {1}{\gamma} [ \mathbf {Q} - (\hat {\mathbf {k}} \cdot \mathbf {Q}) \hat {\mathbf {k}} ] e ^ {[ i \mathbf {k} _ {\parallel} \cdot (\mathbf {r} - \mathbf {r} _ {Q}) + i \gamma | z - z _ {Q} | ]} d \mathbf {k} _ {\parallel}.\tag{10}
$$

In the object space, z is larger than $z _ { Q }$ . Thus, the absolute value is removed and the field can be written as $\mathrm { E q . ~ } ( 5 )$ , with ${ \bf e } ( { \bf k } _ { \parallel } )$ given by 

$$
\mathbf {e} (\mathbf {k} _ {\parallel}) = \frac {i}{8 \pi^ {2}} \frac {1}{\gamma} [ \mathbf {Q} - (\hat {\mathbf {k}} \cdot \mathbf {Q}) \hat {\mathbf {k}} ] \exp [ - i \mathbf {k} \cdot (\mathbf {r} _ {Q} - \mathbf {r} _ {o}) ],\tag{11}
$$

where $\hat { \mathbf { k } } = \mathbf { k } / ( n k _ { 0 } )$ . It is worth noting that this expression is valid no matter the position of the source with respect to the object focal plane. 

Finally, the image field $\mathbf { E } _ { Q , \mathrm { i m } } ( \mathbf { r } _ { \parallel } , z _ { i } )$ is obtained by replacing Eq. (11) in Eq. (3) and in Eq. (7). 

A Matlab code simulating the image of a radiating dipole is provided (see Code 1, Ref. [12]). 

It is worth mentioning that the development of the phase term of the electric field of the dipole in the image space [injecting Eq. (11) in Eq. (6)] enables the derivation of a counterintuitive effect regarding the optimal position of the camera for a defocused object, as demonstrated in Appendix D. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-10-16/79454edb-07dd-45d2-b477-d4dc46f9cffc/d07516e708ff63c2c7e6c4fde8bec06b8737d66d76add15f6ed250d554fcbceb.jpg)



Fig 2 Illustration of a point source Q in a layered medium. One radiated plane wave is transmitted at the interface, $z = 0$ , between two media of refractive indices $n _ { b }$ and n.


## D. Point Source before an Interface

We now consider the common case of a vectorial point source placed in the vicinity ofthe interface oftwo semi-infinite media (see Fig. 2). This configuration typically corresponds to the imaging of a small object in water or air, deposited on a glass coverslip, and imaged by an oil-immersion objective. We introduce $\varepsilon ( z )$ the permittivity of the layered medium before the objective, where $\varepsilon ( z < 0 ) = n _ { b } ^ { 2 }$ and $\varepsilon ( z \ge 0 ) = n ^ { 2 }$ . We recall that the objective is placed somewhere at $z > 0$ . The point source is located at $\mathbf { r } _ { Q }$ with $z _ { Q } < 0$ in the medium ofrefractive index $n _ { b }$ We now define the object space as the domain that lies between $z = 0$ and the objective with refractive index $^ { n , }$ while the image space remains unchanged. As previously, one needs to calculate the plane wave expansion of the field radiated by the point source in the object space, $\mathbf { E } _ { \mathrm { o b } }$ . The field is the solution of the equation 

$$
\nabla \times (\nabla \times \mathbf {E} _ {Q}) - \varepsilon (z) k _ {0} ^ {2} \mathbf {E} _ {Q} = \mathbf {Q} \delta (\mathbf {r} - \mathbf {r} _ {Q}),\tag{12}
$$

which satisfies outgoing boundary condition. The source in the medium of refractive index $n _ { b }$ generates for $0 > z > z _ { Q }$ a sum of plane waves propagating toward positive ${ \mathcal { Z } } ,$ with wavevector ${ \bf k } _ { b } = ( { \bf k } _ { \parallel } , \gamma _ { b } = \sqrt { n _ { b } ^ { 2 } k _ { 0 } ^ { 2 } - k _ { \parallel } ^ { 2 } } )$ and amplitudes ${ \bf e } _ { b } ( { \bf k } _ { \parallel } )$ given by Eq. (11) with n replaced by $n _ { b } .$ Each plane wave is transmitted into the medium of refractive index n at $z = 0$ into a plane wave of amplitude ${ \bf e } ( { \bf k } _ { \parallel } )$ with wavevector $\mathbf { k } = ( \mathbf { k } _ { | | } , \gamma )$ , where $\gamma = \sqrt { n ^ { 2 } k _ { 0 } ^ { 2 } - k _ { \parallel } ^ { 2 } }$ . Introducing the polarization basis $( \pmb { \mathscr { s } } , \pmb { \mathsf { p } } )$ related to k and $( \mathsf { \pmb { s } } _ { b } , \mathsf { \pmb { p } } _ { b } )$ related to $\mathbf { k } _ { b } ,$ , together with the transmission coefficients related to s and p polarizations, $\begin{array} { r } { t _ { s } = \frac { 2 \gamma _ { b } } { \gamma _ { b } + \gamma } , } \end{array}$ $\begin{array} { r }  t _ { \phi } = \frac { 2 n _ { b } n \gamma _ { b } } { n ^ { 2 } \gamma _ { b } + n _ { b } ^ { 2 } \gamma } \ast \ \end{array}$ , respectively, the amplitude of the transmitted plane wave ${ \bf e } ( { \bf k } _ { \parallel } )$ can be written as 

$$
\mathbf {e} (\mathbf {k} _ {\parallel}) = t _ {s} \left(\mathbf {e} _ {b} \left(\mathbf {k} _ {\parallel}\right) \cdot \mathbf {s} _ {b}\right) \mathbf {s} + t _ {p} \left(\mathbf {e} _ {b} \left(\mathbf {k} _ {\parallel}\right) \cdot \mathbf {p} _ {b}\right) \mathbf {p}.\tag{13}
$$

The field in the object space radiated by a source in the vicinity of an interface can thus be written as Eq. (5) with [8,18] 

$$
\mathbf {e} (\mathbf {k} _ {\parallel}) = \frac {i}{8 \pi^ {2}} \frac {1}{\gamma_ {b}} \left(\mathbf {Q} _ {s} + \mathbf {Q} _ {p}\right) \exp \left[ - i \mathbf {k} _ {b} \cdot \mathbf {r} _ {Q} + i \mathbf {k} \cdot \mathbf {r} _ {o} \right],\tag{14}
$$

where 

$$
\mathbf {Q} _ {s} = t _ {s} [ (\mathbf {Q} - (\hat {\mathbf {k}} _ {b} \cdot \mathbf {Q}) \hat {\mathbf {k}} _ {b}) \cdot \mathbf {s} _ {b} ] \mathbf {s},\tag{15}
$$

$$
\mathbf {Q} _ {p} = t _ {p} [ (\mathbf {Q} - (\hat {\mathbf {k}} _ {b} \cdot \mathbf {Q}) \hat {\mathbf {k}} _ {b}) \cdot \mathbf {p} _ {b} ] \mathbf {p}.\tag{16}
$$

To calculate the field at the image plane, it suffices to apply the transformation between e and $\mathbf { e ^ { \prime } }$ given by $\mathrm { E q . ~ } ( 7 )$ 

A Matlab code simulating the computation ofthe image ofa dipole in presence of a substrate is provided (see Code 1, Ref. [12]). 

## 3. ILLUSTRATION WITH A SPHERICAL METAL NANOPARTICLE

To support and illustrate the quantitative model ofthe image of a point source through a microscope, we compare in this section simulated images of nanoparticles to measurements. 

We considered gold nanospheres deposited on a glass sub strate and covered with glycerol using an oil-immersion objective with $M = - 1 0 0$ and $\mathrm { N A } = 1 . { \overset { } { 3 } }$ . Because the refractive indices of oil and glycerol match that of the glass slide, we can assume that the spheres are in a homogeneous medium of refractive index $n = 1 . 5$ . The particles were illuminated by a collimated laser beam at wavelength $\lambda = 5 9 0$ nm, which is assimilated to a linearly x-polarized plane wave propagating along the optical axis. 

We use QLSI [19] to measure the intensity and phase of the image field. More precisely, we measure a complex transmission coefficient, 

$$
t _ {\mathrm{im}} = \sqrt {T (x , y)} e ^ {i \frac {2 \pi}{\lambda} \delta \ell (x, y)},\tag{17}
$$

where $T ( x , y )$ is the intensity ofthe image total field divided by the intensity of the image incident field (i.e., the incident collimated beam). It is usually called the transmittance ofthe object. The length $\delta \ell ( x , y )$ is related to the phase of the field compo nent directed along x and represents the optical path difference (OPD) induced by the object. 

To simulate the transmission coefficient, we introduce $\mathbf { E } _ { \mathrm { e x , o b } } ,$ , the incident field that excites the nanoparticle. In the object space, the latter reads $\mathbf { E } _ { \mathrm { e x , o b } } ( \mathbf { r } _ { \| } , z ) = E _ { 0 }$ exp ink<sub>0</sub>z xˆ. Using Eq. (2), the incident field at the image plane is given by 

$$
\mathbf {E} _ {\mathrm{ex,im}} = \frac {1}{M} \sqrt {n} E _ {0} \exp (- i n k _ {0} z _ {o}) \hat {\mathbf {x}}.\tag{18}
$$

The total field at the image plane is the sum of the incident field and the field scattered by the particle, $\mathbf { E } _ { \mathrm { t o t , i m } } =$ ${ \bf E } _ { \mathrm { e x , i m } } + { \bf E } _ { Q , \mathrm { i m } } .$ The transmission intensity reads $T =$ $E _ { \mathrm { t o t , i m } } ^ { 2 } / E _ { \mathrm { e x , i m } } ^ { 2 } ,$ while the phase is extracted from the complex function, $\dot { \mathbf { E } } _ { \mathrm { t o t , i m } } \cdot \hat { \mathbf { x } } / \mathbf { E } _ { \mathrm { e x , i m } } \cdot \hat { \mathbf { x } }$ 

We now assume that the particle is small enough so that its scattered field is similar to that radiated by a vectorial point source [20], Q, whose amplitude is related to the excitation field at the position of the particle center $\mathbf { r } _ { Q } ,$ with $\mathbf { Q } = k _ { 0 } ^ { 2 } \alpha \mathbf { E } _ { \mathrm { e x , o b } } ( \mathbf { r } _ { \mathbf { Q } } )$ , and α is the complex polarizability of the nanoparticle. In our calculation, the polarizability of the nanoparticle $\begin{array} { r } { \alpha = i \frac { 6 \pi } { n k _ { \alpha } ^ { 3 } } a _ { 1 } \left[ 2 1 \right] } \end{array}$ is determined using Mie theory, where $\scriptstyle a _ { 1 }$ is the dipolar Mie coefficient, given in Appendix C. 

Under these hypotheses, $\mathbf { E } _ { Q , \mathrm { i m } }$ is easily computed using Eqs. (11) and (7). 

Figure 3(a) displays the measured and simulated OPD and intensity images and their profiles for several positions $z _ { Q }$ of the particle with respect to the object focal plane (mechanically displaced by step of about 500 nm). To accurately compare simulations and experiments, three parameters have to be optimized: (i) In theory, for a perfect objective, the phase profile in the microscope objective pupil is supposed to be uniform. The analysis of the experimental images in the Fourier space revealed that the phase profile of the microscope objective pupil was not uniform near the edges of the objective aperture, i.e., at high numerical apertures. This is due to some aberrations created by the optical components and the microscope coverslip, which could be included in the apodization function [15]. Thus, in order to discard these aberrations effects and improve the comparison of our experimental results with simulations, we cropped the field at the microscope entrance pupil calculated from the experimental images at a numerical aperture of 0.9, instead of leaving it at 1.3. (ii) The precision of the measurements of $z _ { Q }$ positions were not enough, as a modification of a few tens of nanometers can already markedly affect the images [22]. For these reasons, the different $z _ { Q }$ have been considered free parameters in the numerical simulations, and the values indicated in Fig. 3 are those that make simulations fit the experimental profiles. (iii) Finally, the nanoparticles we used are supposed to have a diameter of 100 nm, but only in average. We determined the size of the nanoparticle by comparing the experimental and simulated images obtained close to the focus. The best agreement between simulations and theory were obtained for a nanoparticle diameter of 74 nm. Under all these considerations, we have found an excellent quantitative agreement between the experiments and the theory, as shown in Fig. 3. 

In addition, we plot in Fig. 4 the simulation of the scattered intensity along x of a gold nanoparticle as a function of $z = z _ { Q } - z _ { o } ,$ , as well as to the total intensity and OPD images. One can see that the contrast of the OPD image is inverted with respect to a certain axial plane z close to the focus (se dashed line on Fig. 4), while the contrast on the intensity 


(a)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-10-16/79454edb-07dd-45d2-b477-d4dc46f9cffc/24b17dda9bd6ef5a2374f1bbcc76363417e07a126c850a7375fc31215edf84c0.jpg)



Fig 3 (a) OPD and (b) intensity images and the corresponding profiles of a 74 nm gold nanoparticle embedded in a homogeneous medium (n 1.5) at different z positions. The experimental results measured by QLSI are compared to the computed results using the formalism developed in this paper (the homogeneous medium Matlab code provided with this article; see Code 1, Ref. [12]). For both theory and measurements, th illumination wavelength is λ 590 nm, with an objective of 100× magnification and a corrected numerical aperture of 0.9. The scale bar is 2 μm


![image](https://cdn-mineru.openxlab.org.cn/result/2026-10-16/79454edb-07dd-45d2-b477-d4dc46f9cffc/7312993226405bede15b54ffaf3b147330e3cb3ebe4b59f49fe308df962ffe05.jpg)



Fig 4 Simulated images of a single gold nanoparticle (74 nm) in a homogeneous medium $( n = 1 . 5 )$ as a function of $z = z _ { Q } - z _ { o } .$ (a) Normalized magnitude of the scattered field by the nanoparticle. (b) Intensity of the normalized total field $( \mathbf { E } _ { \mathrm { e x , i m } } + \mathbf { E } _ { Q , \mathrm { i m } } )$ and (c) its OPD.


images remains approximately the same. Note that the position ofthe symmetry (or anti-symmetry) plane depends on the characteristics of the nanoparticle (composition, size …) as well as the illumination wavelength (at or off resonance). For example, for a gold nanoparticle illuminated at resonance, the plane is located exactly at the focus plane $z = 0$ 

It should be noted that in the above analysis, we have considered nanoparticles with sizes much smaller than the illumination wavelength in order to satisfy the dipolar approximation. However, our model can be applied to larger objects. In this case, one needs first to calculate the plane-wave expansion of the scattered electromagnetic field in the object space with a solver of Maxwell equations. One can use, for example, the discrete dipole approximation, which consists of discretizing the object into an assembly of interacting dipoles [23]. Then, the sample image is obtained as a coherent sum of the field image ofeach dipole. For a simulator ofmicroscope images based on a rigorous solving of the light–sample interaction, see [24]. 

In summary, we have introduced a model to compute the image, in intensity and phase, of an electric dipole moment through a microscope. The model enables the conception of simple numerical codes, based only on Fourier transforms, that do not require the calculation of integrals. Note that Eq. (6) even enables the calculation of the image of any object, not only dipoles, from its decomposition in plane waves. A good quantitative agreement was obtained with phase and intensity measurements of a light beam obtained with QLSI. Our study paves the way for the quantitative analysis of optical properties of small scatterers, such as molecules or nanoparticles, and in general should be useful for metrology applications involving microscopy imaging. 

Open-source code that permits calculation of the image of a dipole through a microscope is provided (see Code 1, Ref. [12]). 

## APPENDIX A: ENERGY CONSERVATION BEFORE AND AFTER THE IMAGING SYSTEM

We demonstrate here how the factor $\begin{array} { r } { A = \frac { 1 } { M } \sqrt { \frac { \gamma } { \gamma ^ { \prime } } } } \end{array}$ in Eqs. (2) and (6) is obtained based on energy conservation considerations. The field in the object space is given by 

$$
\mathbf {E} _ {\mathrm{ob}} (\mathbf {r}) = \int \mathbf {e} (\mathbf {k} _ {\parallel}) \exp [ i \mathbf {k} _ {\parallel} \cdot \mathbf {r} _ {\parallel} + i \gamma z ] \mathrm{d} \mathbf {k} _ {\parallel},\tag{A1}
$$

while the field at the image space is given by 

$$
\mathbf {E} _ {\mathrm{im}} (\mathbf {r}) = \int A \mathbf {e} ^ {\prime} (\mathbf {k} _ {\parallel}) \exp [ i \mathbf {k} _ {\parallel} ^ {\prime} \cdot \mathbf {r} _ {\parallel} + i \gamma^ {\prime} z ] \mathrm{d} \mathbf {k} _ {\parallel},\tag{A2}
$$

with $| \mathbf { e } ( \mathbf { k } _ { \| } ) | = | \mathbf { e } ^ { \prime } ( \mathbf { k } _ { \| } ) |$ . A is a factor that ensures the energy conservation between the object space and image space. To determine this factor, we calculate the transmitted power through an infinite transverse plane S in the object space and image space: 

$$
P _ {s} = \frac {1}{2} \int \operatorname{Re} [ \mathbf {E} \times \mathbf {H} ^ {*} ] \cdot \hat {\mathbf {z}} d \mathbf {r} _ {\parallel},\tag{A3}
$$

where H, the magnetic field, is expressed as a function of th electric field as $\mathbf { H } = \nabla \times \mathbf { E } / i \mu _ { 0 } \omega$ 

Let us first calculate the transmitted power in the object space. The magnetic field is 

$$
\mathbf {H} _ {\mathrm{ob}} (\mathbf {r}) = \frac {1}{i \mu_ {0} \omega} \int i \mathbf {k} \times \mathbf {e} (\mathbf {k} _ {\parallel}) \exp [ i \mathbf {k} _ {\parallel} \cdot \mathbf {r} _ {\parallel} + i \gamma z ] \mathrm{d} \mathbf {k} _ {\parallel},\tag{A4}
$$

where $\mu _ { 0 }$ is the vacuum permeability and $\omega = 2 \pi / \lambda$ . Thus, the transmitted power at the object plane reads, with obvious notations, 

$$
\begin{array}{l} P _ {\mathrm{sob}} = \frac {1}{2 \mu_ {0} \omega} \int \operatorname{Re} \left[ \int \mathbf {e} ^ {*} (\mathbf {k} _ {\parallel 1}) \exp [ - i \mathbf {k} _ {\parallel 1} \cdot \mathbf {r} _ {\parallel} - i \gamma_ {1} z ] \mathrm{d} \mathbf {k} _ {\parallel 1} \right. \\ \times \int \mathbf {k} _ {2} \times \mathbf {e} (\mathbf {k} _ {\parallel 2}) \exp [ i \mathbf {k} _ {\parallel 2} \cdot \mathbf {r} _ {\parallel} + i \gamma_ {2} z ] \mathrm{d} \mathbf {k} _ {\parallel 2} ] \mathrm{d} \mathbf {r} _ {\parallel} \cdot \hat {\mathbf {z}}. \end{array}\tag{A5}
$$

Exchanging the order of the spatial and frequency integrations, one obtains 

$$
\begin{array}{l} P _ {\mathrm{sob}} = \frac {1}{2 \mu_ {0} \omega} \mathrm{Re} \bigg [ \iint \mathbf {e} ^ {*} (\mathbf {k} _ {\parallel 1}) \times (\mathbf {k} _ {2} \times \mathbf {e} (\mathbf {k} _ {\parallel 2})) \\ \quad \times \exp [ i (\gamma_ {2} - \gamma_ {1}) z ] \\ \quad \times \left[ \int \exp [ i (\mathbf {k} _ {\parallel 2} - \mathbf {k} _ {\parallel 1}) \cdot \mathbf {r} _ {\parallel} ] \mathrm{d} \mathbf {r} _ {\parallel} \cdot \hat {\mathbf {z}} \right] \mathrm{d} \mathbf {k} _ {\parallel 1} \mathrm{d} \mathbf {k} _ {\parallel 2} \bigg ]. \end{array}\tag{A6}
$$

Integrating the spatial variables, we get 

$$
\begin{array}{l} P _ {\text { sob }} = \frac {2 \pi^ {2}}{\mu_ {0} \omega} \mathrm{Re} \bigg [ \iint \mathbf {e} ^ {*} (\mathbf {k} _ {\parallel 1}) \times (\mathbf {k} _ {2} \times \mathbf {e} (\mathbf {k} _ {\parallel 2})) \\ \cdot \hat {\mathbf {z}} \exp [ i (\gamma_ {2} - \gamma_ {1}) z ] \delta (\mathbf {k} _ {\parallel 2} - \mathbf {k} _ {\parallel 1}) d \mathbf {k} _ {\parallel 1} d \mathbf {k} _ {\parallel 2} \bigg ], \end{array}\tag{A7}
$$

where δ is the Dirac function. Knowing that ${ \bf e } ^ { * } ( { \bf k } _ { | | } ) \times ( { \bf k } \times$ ${ \bf e } ( { \bf k } _ { \parallel } ) ) = | { \bf e } ( { \bf k } _ { \parallel } ) | ^ { 2 } { \bf k }$ k and k $\hat { \mathbf { z } } = \boldsymbol { \gamma }$ , the transmitted power at the object plane reads 

$$
P _ {\mathrm{sob}} = \frac {2 \pi^ {2}}{\mu_ {0} \omega} \int | \mathbf {e} (\mathbf {k} _ {\parallel}) | ^ {2} \gamma \mathrm{d} \mathbf {k} _ {\parallel}.\tag{A8}
$$

In the same manner, we calculate the transmitted power at the image plane: 

$$
\begin{array}{l} P _ {\text { sim }} = \frac {1}{2 \mu_ {0} \omega} \operatorname{Re} \left[ \iint A ^ {2} \mathbf {e} ^ {\prime *} (\mathbf {k} _ {\| 1}) \right. \\ \times \left(\mathbf {k} _ {2} ^ {\prime} \times \mathbf {e} ^ {\prime} (\mathbf {k} _ {\| 2})\right) \left[ \int \exp \left[ i \left(\frac {\mathbf {k} _ {\| 2} - \mathbf {k} _ {\| 1}}{M}\right) \cdot \mathbf {r} _ {\|} \right] d \mathbf {r} _ {\|}. \hat {\mathbf {z}} \right] \\ \times \exp [ i (\gamma_ {2} - \gamma_ {1}) z ] d \mathbf {k} _ {\| 1} d \mathbf {k} _ {\| 2}. \end{array} \tag {A}\tag{A9}
$$

The integration over the spatial coordinates is performed using the change of variables, $\dot { \mathbf { u } } = \mathbf { r } _ { \| } / M , ~ \mathrm { d } \mathbf { r } _ { \| } = \dot { M ^ { 2 } } \mathrm { d } \mathbf { u }$ , Eq. (A9): 

$$
\begin{array}{l} P _ {\text { sim }} = \frac {2 \pi^ {2}}{\mu_ {0} \omega} \operatorname{Re} \left[ \iint A ^ {2} \mathbf {e} ^ {\prime *} (\mathbf {k} _ {\| 1}) \times (\mathbf {k} _ {2} ^ {\prime} \times \mathbf {e} ^ {\prime} (\mathbf {k} _ {\| 2})) \right. \\ \cdot \hat {\mathbf {z}} \exp [ i (\gamma_ {2} - \gamma_ {1}) z ] M ^ {2} \delta (\mathbf {k} _ {\| 2} - \mathbf {k} _ {\| 1}) d \mathbf {k} _ {\| 1} d \mathbf {k} _ {\| 2} \Bigg ]. \end{array}\tag{A10}
$$

The simplification of Eq. (A10) yields 

$$
P _ {\mathrm{sim}} = \frac {2 \pi^ {2}}{\mu_ {0} \omega} M ^ {2} \int A ^ {2} | \mathbf {e} (\mathbf {k} _ {\parallel}) | ^ {2} \gamma^ {\prime} \mathrm{d} \mathbf {k} _ {\parallel}.\tag{A11}
$$

The energy conservation between the object space and image space states that $P _ { \mathrm { s o b } } = P _ { \mathrm { s i m } } ;$ thus, one can extract the factor A as 

$$
A = \frac {1}{M} \sqrt {\frac {\gamma}{\gamma^ {\prime}}}.\tag{A12}
$$

## APPENDIX B: ROTATION MATRIX

The rotation matrix $\mathbf { R } ( \mathbf { k } _ { | | } )$ is given by 

$$
\begin{array}{l} \mathbf {R} (\mathbf {k} _ {\parallel}) \\ = \left( \begin{array}{c c c} u ^ {2} (1 - \cos \theta) + \cos \theta & u u _ {y} (1 - \cos \theta) & u _ {y} \sin \theta \\ u u _ {y} (1 - \cos \theta) & u _ {y} ^ {2} (1 - \cos \theta) + \cos \theta & - u \sin \theta \\ - u _ {y} \sin \theta & u \sin \theta & \cos \theta \end{array} \right), \end{array}\tag{B1}
$$

where $\begin{array} { r } { \mathbf { u } = \frac { \hat { \mathbf { k } } \times \mathbf { z } } { | \hat { \mathbf { k } } \times \mathbf { z } | } } \end{array}$ is the rotation axis. Note that u has no component along the z direction. θ is defined as cos $\theta = \hat { \mathbf { k } } \cdot \hat { \mathbf { k } } ^ { \prime }$ and sin $\theta = | \hat { \mathbf { k } } \times \hat { \mathbf { k } } ^ { \prime } |$ as shown in Fig. 5. 

## APPENDIX C: DIPOLAR MIE COEFFICIENT

The dipolar Mie coefficient for a nanosphere of radius r having refractive index $n _ { s }$ embedded in a medium with refractive index n reads 

$$
a _ {1} = \frac {m \psi_ {1} (w) \psi_ {1} ^ {\prime} (v) - \psi_ {1} (v) \psi_ {1} ^ {\prime} (w)}{m \psi_ {1} (w) \zeta_ {1} ^ {\prime} (v) - \zeta_ {1} (v) \psi_ {1} ^ {\prime} (w)},\tag{C1}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-10-16/79454edb-07dd-45d2-b477-d4dc46f9cffc/e175e3cbd7dce7b29ac29f6d7227930ff79ee4161a8c9d8728b62571f6467741.jpg)



Fig 5 Illustration ofthe rotation axis and rotation angle used in the rotation matrix $\mathbf { R } ( \mathbf { k } _ { | | } )$ .


where $m = n _ { s } / n , \ \nu = k r ,$ and $w = m \nu$ . The functions $\psi _ { 1 } , \zeta$ and their derivatives $\psi _ { 1 } ^ { \prime } , \zeta _ { 1 } ^ { \prime }$ are given by 

$$
\psi_ {1} (x) = \sin (x) / x - \cos (x),\tag{C2}
$$

$$
\zeta_ {1} (x) = \sin (x) / x - i (\cos (x) / x + \sin (x)),\tag{C3}
$$

$$
\psi_ {1} ^ {\prime} (x) = \sin (x) - \psi_ {1} (x) / x,\tag{C4}
$$

$$
\zeta_ {1} ^ {\prime} (x) = \sin (x) - i \cos (x) - \zeta_ {1} (x) / x.\tag{C5}
$$

## APPENDIX D: DEFOCUSING THE OBJECTIVE OR THE CAMERA: IS IT THE SAME?

When a point source is located at the object focal plane, one knows that the image of the source (i.e., the smallest intensity spot) will be located at the image focal plane. In this appendix, we discuss the case of a source placed before or after the object focal plane. Is it possible to displace the detector out of the image focal plane to retrieve a focused spot? In other words, can we find a position where all the plane waves are in phase in the image space? To this aim, we develop the phase term ofthe plane waves emitted by a point source that are propagating in the image space. The phase obtained by injecting Eq. (11) into Eq. (6) i 

$$
\begin{array}{l} \exp [ i \mathbf {k} ^ {\prime} \cdot (\mathbf {r} - \mathbf {r} _ {i}) - i \mathbf {k} \cdot (\mathbf {r} _ {Q} - \mathbf {r} _ {o}) ] \\ = \exp \bigg [ i \mathbf {k} _ {\parallel} \cdot \bigg (\frac {\mathbf {r} _ {\parallel}}{M} - \mathbf {r} _ {\parallel , Q} \bigg) + i \gamma^ {\prime} (z - z _ {i}) - i \gamma (z _ {Q} - z _ {o}) \bigg ]. \end{array}\tag{D1}
$$

Equation (D1) shows that if the point source is placed at the object focal plane, $z _ { Q } = z _ { o } ,$ , then all the plane waves in the image space converge (i.e., have the same phase) at $z = z _ { i }$ and $\mathbf { r } _ { \parallel } = M \mathbf { r } _ { \parallel , Q } .$ . As expected, we obtain a focused spot at the image focal plane, and the transverse position has been magnified by M. 

On the contrary, if the point source is placed before or after the object plane, $z _ { Q } \neq z _ { o } ,$ , then there is no position in the image space where all the plane waves interfere constructively. Thus, if the point source is out of the object focal plane, one cannot retrieve a focused spot in the image space that is as small as that obtained when the point source is in the object focal plane. 

Yet, in the case of the paraxial approximation where $k _ { \parallel } \ll n k _ { 0 }$ , using $\prime \approx n k _ { 0 } - k _ { \parallel } ^ { 2 } / ( \hat { 2 n } k _ { 0 } )$ and $\dot { \gamma ^ { \prime } } \approx k _ { 0 } - k _ { \parallel } ^ { 2 } / ( 2 M ^ { 2 } k _ { 0 } )$ the phase term ofthe plane waves $\left[ \mathrm { E q . } \left( \mathrm { D 1 } \right) \right]$ can be written, to a constant, as 

$$
\exp \left[ i \mathbf {k} _ {\parallel} \cdot \left(\frac {\mathbf {r} _ {\parallel}}{M} - \mathbf {r} _ {\parallel , Q}\right) - i \frac {k _ {\parallel} ^ {2}}{2 k _ {0}} \left[ \frac {z - z _ {i}}{M ^ {2}} - \frac {z _ {Q} - z _ {o}}{n} \right] \right].\tag{D2}
$$

Under this approximation, all the plane waves interfere constructively when $\hat { \mathbf { r } } = M \mathbf { r } _ { \parallel , Q } + [ z _ { i } + \hat { M } ^ { 2 } ( z _ { Q } - z _ { o } ) / n ] \hat { \mathbf { z } }$ . Thus, under paraxial approximation, if a point source is located at a distance $d$ from the object focal plane, one retrieves its image $( \mathrm { i . e . , }$ the best possible focused spot) at a distance $M ^ { 2 } / n d$ from the image focal plane, the magnification in the transverse direction remaining M. 

Funding. Agence Nationale de la Recherche (ANR) (ANR-16-CE09-0002). 

Acknowledgment. This work was supported by the French Agence Nationale de la Recherche (ANR NANOHEATERS, Grant ANR-16-CE09-0002). 

## REFERENCES



1. F. A. Merchant, “Three-dimensional imaging,” in Microscope Image Processing, Q. Wu, F. A. Merchant, and K. R. Castleman, eds. (Academic, 2008), Chap. 14, pp. 329–399. 





2. E. Dusch, T. Dorval, N. Vincent, M. Wachsmuth, and A. Genovesio, “Three-dimensional point spread function model for line-scanning confocal microscope with high-aperture objective,” J. Microsc. 228, 132–138 (2007). 





3. B. Richards and E. Wolf, “Electromagnetic diffraction in optica systems. II. Structure of the image field in an aplanatic system,” Proc. R. Soc. London Ser. A 253, 358–379 (1959). 





4. N. Streibl, “Three-dimensional imaging by a microscope,” J. Opt. Soc. Am. A 2, 121–127 (1985). 





5. F. Aguet, D. V. D. Ville, and M. Unser, “Model-based 2.5-D deconvolution for extended depth of field in brightfield microscopy,” IEEE Trans. Image Process. 17, 1144–1153 (2008). 





6. S. Stallinga and B. Rieger, “Accuracy of the Gaussian point spread function model in 2D localization microscopy,” Opt. Express 18, 24461–24476 (2010). 





7. J. Enferlein, “Theoretical study of detection of a dipole emitter through an objective with high numerical aperture,” Opt. Lett. 25, 634–636 (2000). 





8. L. Novotny and B. Hecht, Principle of Nano-Optics (Cambridge University, 2006). 





9. A. S. Backer and W. E. Moerner, “Extending single-molecule microscopy using optical Fourier processing,” J. Phys. Chem. B 118, 8313–8329 (2014). 





10. T. Zhang, Y. Ruan, G. Maire, D. Sentenac, A. Talneau, K. Belkebir, P. C. Caumet, and A. Sentenac, “Full-polarized tomographic diffraction microscopy achieves a resolution about one-fourth of the wavelength,” Phys. Rev. Lett. 111, 243904 (2013). 





11. B. J. Davis, P. S. Carney, and R. Bhargava, “Theory of midinfrared absorption microspectroscopy: I. Homogeneous samples,” Anal. Chem. 82, 3474–3486 (2010). 





12. S. Khadir, P. C. Chaumet, G. Baffou, and A. Sentenac, “Open-source code that permits calculation of the image of a dipole through a 





microscope,” figshare (2019), https://doi.org/10.6084/m9.figshare .7704494. 





13. E. Abbe, “VII.—On the estimation of aperture in the microscope,” J. R. Microsc. Soc. 1, 388–423 (1981) 





14. C. J. R. Sheppard and K. G. Larkin, “Effect of numerical aperture on interference fringe spacing,” Appl. Opt. 34, 4731–4734 (1995). 





15. H. Zhou, M. Gu, and C. J. R. Sheppard, “Investigation of aberration measurement in confocal microscopy,” J. Modern Opt. 42, 627–638 (1995). 





16. J. A. Kurvits, M. Jiang, and R. Zia, “Comparative analysis of imaging configurations and objectives for Fourier microscopy,” J. Opt. Soc. Am. A 32, 2082–2092 (2015) 





17. L. Tsang, J. A. Kong, K. H. Ding, and C. O. Ao, Scattering of Electromagnetic Waves (Wiley-Interscience, 2001). 





18. G. S. Agarwal, “Quantum electrodynamics in the presence of dielectrics and conductors. I Electromagnetic-field response functions and black-body fluctuations in finite geometry,” Phys. Rev. A 11, 230–242 (1975). 





19. S. Khadir, P. Bon, D. Vignaud, E. Galopin, N. McEvoy, D. McCloskey, S. Monneret, and G. Baffou, “Optical imaging and characterization of graphene and other 2D materials using quantitative phase microscopy,” ACS Photon. 4, 3130–3139 (2017). 





20. P. C. Chaumet, A. Rahmani, F. de Fornel, and J.-P. Dufour, “Evanescent light scattering: the validity of the dipole approximation, Phys. Rev. B 58, 2310–2315 (1998). 





21. R. Colom, A. Devilez, N. Bonod, and B. Stout, “Optimal interactions of light with magnetic and electric resonant particles,” Phys. Rev. B 93, 045427 (2016). 





22. P. Bon, N. Bourg, S. Lécart, S. Monneret, E. Fort, J. Wenger, and S. Lévêque-Fort, “Three-dimensional nanometre localization of nanopar ticles to enhance super-resolution microscopy,” Nat. Commun. 6, 7764 (2015). 





23. P. C. Chaumet, A. Sentenac, and A. Rahmani, “Coupled dipole method for scatterers with large permittivity,” Phys. Rev. E 70, 036606 (2004). 





24. https://github.com/pchaumet/IF-DDA/. 

