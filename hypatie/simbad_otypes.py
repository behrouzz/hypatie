text = """00.00.00.00,?,Object of unknown nature
00.02.00.00,ev,transient event
01.00.00.00,Rad,Radio-source
01.02.00.00,mR,metric Radio-source
01.04.00.00,cm,centimetric Radio-source
01.06.00.00,mm,millimetric Radio-source
01.08.00.00,smm,sub-millimetric source
01.11.00.00,HI,HI (21cm) source
01.12.00.00,rB,radio Burst
01.14.00.00,Mas,Maser
02.00.00.00,IR,Infra-Red source
02.02.00.00,FIR,Far-Infrared source
02.03.00.00,MIR,Mid-Infrared source
02.04.00.00,NIR,Near-Infrared source
03.00.00.00,red,Very red source
03.03.00.00,ERO,Extremely Red Object
04.00.00.00,blu,Blue object
05.00.00.00,UV,UV-emission source
06.00.00.00,X,X-ray source
06.02.00.00,UX?,Ultra-luminous X-ray candidate
06.10.00.00,ULX,Ultra-luminous X-ray source
07.00.00.00,gam,gamma-ray source
07.03.00.00,gB,gamma-ray Burst
08.00.00.00,err,"Not an object (error, artefact, ...)"
09.00.00.00,grv,Gravitational Source
09.03.00.00,Lev,(Micro)Lensing Event
09.06.00.00,LS?,Possible gravitational lens System
09.07.00.00,Le?,Possible gravitational lens
09.08.00.00,LI?,Possible gravitationally lensed image
09.09.00.00,gLe,Gravitational Lens
09.11.00.00,gLS,Gravitational Lens System (lens+images)
09.12.00.00,GWE,Gravitational Wave Event
10.00.00.00,..?,Candidate objects
10.01.00.00,G?,Possible Galaxy
10.02.00.00,SC?,Possible Supercluster of Galaxies
10.03.00.00,C?G,Possible Cluster of Galaxies
10.04.00.00,Gr?,Possible Group of Galaxies
10.06.00.00,As?,
10.11.00.00,**?,Physical Binary Candidate
10.11.01.00,EB?,Eclipsing Binary Candidate
10.11.10.00,Sy?,Symbiotic Star Candidate
10.11.11.00,CV?,Cataclysmic Binary Candidate
10.11.11.06,No?,Nova Candidate
10.11.12.00,XB?,X-ray binary Candidate
10.11.12.02,LX?,Low-Mass X-ray binary Candidate
10.11.12.03,HX?,High-Mass X-ray binary Candidate
10.12.00.00,Pec?,Possible Peculiar Star
10.12.01.00,Y*?,Young Stellar Object Candidate
10.12.02.00,pr?,Pre-main sequence Star Candidate
10.12.02.03,TT?,T Tau star Candidate
10.12.03.00,C*?,Possible Carbon Star
10.12.04.00,S*?,Possible S Star
10.12.05.00,OH?,Possible Star with envelope of OH/IR type
10.12.07.00,WR?,Possible Wolf-Rayet Star
10.12.08.00,Be?,Possible Be Star
10.12.09.00,Ae?,Possible Herbig Ae/Be Star
10.12.11.00,HB?,Possible Horizontal Branch Star
10.12.11.02,RR?,Possible Star of RR Lyr type
10.12.11.03,Ce?,Possible Cepheid
10.12.11.07,WV?,Possible Variable Star of W Vir type
10.12.12.00,RB?,Possible Red Giant Branch star
10.12.13.00,sg?,Possible Supergiant star
10.12.13.03,s?r,Possible Red supergiant star
10.12.13.04,s?y,Possible Yellow supergiant star
10.12.13.05,s?b,Possible Blue supergiant star
10.12.14.00,AB?,Asymptotic Giant Branch Star candidate
10.12.14.01,LP?,Long Period Variable candidate
10.12.14.02,Mi?,Mira candidate
10.12.15.00,pA?,Post-AGB Star Candidate
10.12.16.00,BS?,Candidate blue Straggler Star
10.12.17.00,HS?,Hot subdwarf candidate
10.12.18.00,WD?,White Dwarf Candidate
10.12.20.00,N*?,Neutron Star Candidate
10.12.22.00,BH?,Black Hole Candidate
10.12.23.00,SN?,SuperNova Candidate
10.12.24.00,LM?,Low-mass star candidate
10.12.26.00,BD?,Brown Dwarf Candidate
12.00.00.00,mul,Composite object
12.01.00.00,reg,Region defined in the sky
12.01.05.00,vid,Underdense region of the Universe
12.02.00.00,SCG,Supercluster of Galaxies
12.03.00.00,ClG,Cluster of Galaxies
12.04.00.00,GrG,Group of Galaxies
12.04.05.00,CGG,Compact Group of Galaxies
12.05.00.00,PaG,Pair of Galaxies
12.05.05.00,IG,Interacting Galaxies
12.09.00.00,C?*,Possible (open) star cluster
12.10.00.00,Gl?,Possible Globular Cluster
12.11.00.00,Cl*,Cluster of Stars
12.11.01.00,GlC,Globular Cluster
12.11.02.00,OpC,Open (galactic) Cluster
12.12.00.00,As*,Association of Stars
12.12.01.00,St*,Stellar Stream
12.12.02.00,MGr,Moving Group
12.13.00.00,**,Double or multiple star
12.13.01.00,EB*,Eclipsing binary
12.13.01.01,Al*,Eclipsing binary of Algol type
12.13.01.02,bL*,Eclipsing binary of beta Lyr type
12.13.01.03,WU*,Eclipsing binary of W UMa type
12.13.02.00,SB*,Spectroscopic binary
12.13.05.00,El*,Ellipsoidal variable Star
12.13.10.00,Sy*,Symbiotic Star
12.13.11.00,CV*,Cataclysmic Variable Star
12.13.11.02,DQ*,CV DQ Her type (intermediate polar)
12.13.11.03,AM*,CV of AM Her type (polar)
12.13.11.05,NL*,Nova-like Star
12.13.11.06,No*,Nova
12.13.11.07,DN*,Dwarf Nova
12.13.12.00,XB*,X-ray Binary
12.13.12.02,LXB,Low Mass X-ray Binary
12.13.12.03,HXB,High Mass X-ray Binary
13.00.00.00,ISM,Interstellar matter
13.01.00.00,PoC,Part of Cloud
13.02.00.00,PN?,Possible Planetary Nebula
13.03.00.00,CGb,Cometary Globule
13.04.00.00,bub,Bubble
13.06.00.00,EmO,Emission Object
13.08.00.00,Cld,Cloud
13.08.03.00,GNe,Galactic Nebula
13.08.04.00,BNe,Bright Nebula
13.08.06.00,DNe,Dark Cloud (nebula)
13.08.07.00,RNe,Reflection Nebula
13.08.12.00,MoC,Molecular Cloud
13.08.12.03,glb,Globule (low-mass dark cloud)
13.08.12.06,cor,Dense core
13.08.12.08,SFR,Star forming region
13.08.13.00,HVC,High-velocity Cloud
13.09.00.00,HII,HII (ionized) region
13.10.00.00,PN,Planetary Nebula
13.11.00.00,sh,HI shell
13.12.00.00,SR?,SuperNova Remnant Candidate
13.13.00.00,SNR,SuperNova Remnant
13.14.00.00,cir,CircumStellar matter
13.14.01.00,of?,Outflow candidate
13.14.15.00,out,Outflow
13.14.16.00,HH,Herbig-Haro Object
14.00.00.00,*,Star
14.01.00.00,*iC,Star in Cluster
14.02.00.00,*iN,Star in Nebula
14.03.00.00,*iA,Star in Association
14.05.00.00,V*?,Star suspected of Variability
14.06.00.00,Pe*,Peculiar Star
14.06.01.00,HB*,Horizontal Branch Star
14.06.02.00,Y*O,Young Stellar Object
14.06.02.04,Ae*,Herbig Ae/Be star
14.06.05.00,Em*,Emission-line Star
14.06.05.03,Be*,Be Star
14.06.06.00,BS*,Blue Straggler Star
14.06.10.00,RG*,Red Giant Branch star
14.06.12.00,AB*,Asymptotic Giant Branch Star (He-burning)
14.06.12.03,C*,Carbon Star
14.06.12.06,S*,S Star
14.06.13.00,sg*,Evolved supergiant star
14.06.13.03,s*r,Red supergiant star
14.06.13.04,s*y,Yellow supergiant star
14.06.13.05,s*b,Blue supergiant star
14.06.14.00,HS*,Hot subdwarf
14.06.15.00,pA*,Post-AGB Star (proto-PN)
14.06.16.00,WD*,White Dwarf
14.06.17.00,LM*,Low-mass star (M<1solMass)
14.06.18.00,BD*,Brown Dwarf (M<0.08solMass)
14.06.19.00,N*,Confirmed Neutron Star
14.06.23.00,OH*,OH/IR star
14.06.25.00,pr*,Pre-main sequence Star
14.06.25.03,TT*,T Tau-type Star
14.06.30.00,WR*,Wolf-Rayet Star
14.07.00.00,PM*,High proper-motion Star
14.08.00.00,HV*,High-velocity Star
14.09.00.00,V*,Variable Star
14.09.01.00,Ir*,Variable Star of irregular type
14.09.01.01,Or*,Variable Star of Orion Type
14.09.03.00,Er*,Eruptive variable Star
14.09.03.04,RC*,Variable Star of R CrB type
14.09.03.05,RC?,Variable Star of R CrB type candiate
14.09.04.00,Ro*,Rotationally variable Star
14.09.04.01,a2*,Variable Star of alpha2 CVn type
14.09.04.03,Psr,Pulsar
14.09.04.04,BY*,Variable of BY Dra type
14.09.04.05,RS*,Variable of RS CVn type
14.09.05.00,Pu*,Pulsating variable Star
14.09.05.02,RR*,Variable Star of RR Lyr type
14.09.05.03,Ce*,Cepheid variable Star
14.09.05.05,dS*,Variable Star of delta Sct type
14.09.05.06,RV*,Variable Star of RV Tau type
14.09.05.07,WV*,Variable Star of W Vir type
14.09.05.08,bC*,Variable Star of beta Cep type
14.09.05.09,cC*,Classical Cepheid (delta Cep type)
14.09.05.10,gD*,Variable Star of gamma Dor type
14.09.05.11,SX*,Variable Star of SX Phe type (subdwarf)
14.09.06.00,LP*,Long-period variable star
14.09.06.01,Mi*,Variable Star of Mira Cet type
14.09.08.00,SN*,SuperNova
14.14.00.00,su*,Sub-stellar object
14.14.02.00,Pl?,Extra-solar Planet Candidate
14.14.10.00,Pl,Extra-solar Confirmed Planet
15.00.00.00,G,Galaxy
15.01.00.00,PoG,Part of a Galaxy
15.02.00.00,GiC,Galaxy in Cluster of Galaxies
15.02.02.00,BiC,Brightest galaxy in a Cluster (BCG)
15.03.00.00,GiG,Galaxy in Group of Galaxies
15.04.00.00,GiP,Galaxy in Pair of Galaxies
15.05.00.00,HzG,Galaxy with high redshift
15.07.00.00,rG,Radio Galaxy
15.08.00.00,H2G,HII Galaxy
15.09.00.00,LSB,Low Surface Brightness Galaxy
15.10.00.00,AG?,Possible Active Galaxy Nucleus
15.10.07.00,Q?,Possible Quasar
15.10.11.00,Bz?,Possible Blazar
15.10.17.00,BL?,Possible BL Lac
15.11.00.00,EmG,Emission-line galaxy
15.12.00.00,SBG,Starburst Galaxy
15.13.00.00,bCG,Blue compact Galaxy
15.14.00.00,LeI,Gravitationally Lensed Image
15.14.01.00,LeG,Gravitationally Lensed Image of a Galaxy
15.14.07.00,LeQ,Gravitationally Lensed Image of a Quasar
15.15.00.00,AGN,Active Galaxy Nucleus
15.15.01.00,LIN,LINER-type Active Galaxy Nucleus
15.15.02.00,SyG,Seyfert Galaxy
15.15.02.01,Sy1,Seyfert 1 Galaxy
15.15.02.02,Sy2,Seyfert 2 Galaxy
15.15.03.00,Bla,Blazar
15.15.03.01,BLL,BL Lac - type object
15.15.03.02,OVV,Optically Violently Variable object
15.15.04.00,QSO,Quasar"""
