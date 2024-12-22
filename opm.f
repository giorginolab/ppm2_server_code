c                *** opm.f ***
c     Universal version with  new solvation model (transmemberane or
c     peripheral proteins and small molecules)
c
c     Modules: solva.f, min.f, rmsd.f, deftm.f, hbcor.f, profile.f,
c     read_small.f, watface.f, find_segm.f, locate.f, orient.f, readpdb.f,
c     tilting.f
c
c                       Method comments
c                       ---------------
c         The initial orientation is automatically determined by the method of
c         Tusnady. This is done with modified old energy functions (ionizable
c         groups are treated using ASA-based model) in 'orient' module with 
c         recangular rather than sigmoidal boundary. Then, the initial
c         TM secondary stuctures are automatically calculated.
c
c         The "small molecule" option is defined by record "Small molecule" in
c         PDB file. Small molecules should be prepared with hydrogens 
c         and additional information about dipole moments and pKa assigned to
c         specific groups (see comments in read_small.f). They currently
c         include only C,N,O,S atoms.
c
c                     Treatment of heteroatoms
c                     -------------------------
c      1. Solvent molecules ('namsolv' list) are always excluded from 
c         the calculations but re-included in the final coordinate file.
c      2. Cofactors (everything not in residue library) are excluded
c         if ihetero=0 but re-included in final coordinate file. 
c         List of residues 'namwork' serves only to define what is treated
c         heteroi-groups in output PDB file (out.pdb)
c      3. Cofactors are treated as amino acid residues (based on their 
c         residue numbers. Heteroatoms with C,O,N,S in 1st position are 
c         treated as Csp2, OH, NH, and S, respectively. Others are
c         assigned as OH. No electrostatic contributions for charged
c         heteroatoms. Cofactors must be included in the library of
c         residues. More atom types should be added (esp. for ions).
c
c                 Conventions and approximations
c                 ------------------------------
c      1. Some residues in PDB files are incomplete, which leads
c         to omitting contributions of their polar atoms, dipoles, and 
c         ionizable groups. Ideally, they must be reconstructed.
c      2. Electrostatic contributions for dipoles and charged groups are 
c         equally divided between participating atoms. Dipole moments
c         are taken for unionized states. Choice of ionized or unionized
c         state (whichever was of lower energy), rtaher than considering
c         the ionized-unionized equiliblium. No smoothening for ionization
c         energy (eioniz).
c      3. H-bonds are used only for OH groups of S,T and Y. They are 
c         identified as in Quanta (<3.5 A distance and bond angles >90 degrees).
c      4. Topology is defined for N-nerminus of 1st subunit in PDB file
c
      character*80 pdbtempl,line,title,ssfile,ssout,pdbout
      character*3 itopo
c
      open (23,file='datasub1')
      open (24,file='datapar1')
c
c     ihetero= 0 - all hetero-groups are excluded
c              1 - all hetero-groups are included except ones
c                  defined as solvent
c
c     model  = 1 - peripheral protein with fixed membrane thickness of 30 A
c              2 - TM protein with variable membrane thickness
c                  (hydrophobic  matching condition depends on type of
c                  sec. structure)
c              3 - peripheral protein (no secondary structure) with variable
c                  membrane thickness (matching condition for alpha-hel. protein)
c
c     ifunc  = 1 - old energy functions with simplified (ASA-dependent) 
c                  energy for ionizable groups and sigmoidal polarity profile
c                  Old functions are used for initial calculation of transmembrane
c                  domain and transfer energy contributions Rfor individual secondary
c                  structures
c              2 - new energy functions
c
      do ipdb=1,150000
         iprint=1
         itopo='in '
c
c     Read names of PDB and .cor files:
c
         read (*,'(i2,1x,a3,1x,a80)',end=40) 
     *     ihetero,itopo,pdbtempl
         write (*,'(a80)') pdbtempl
c
         j00=index(pdbtempl,'.pdb')-1
         if(j00.le.2) then
            write (*,'(''wrong name of PDB file:'',a50)') pdbtempl
            stop
         end if
         pdbout=pdbtempl(1:j00)//'out.pdb'
c
         if(ihetero.eq.1) then
            write (*,'(''heteroatoms included'')')
         else
            write (*,'(''heteroatoms excluded'')')
         end if
c
c        2. Calculate position of the protein
c      
         call immers (pdbtempl,iprint,ihetero,itopo,pdbout)
c
      end do
   40 continue
      stop
      end
c
      subroutine immers (pdbtempl,iprint,ihetero,itopo,pdbout)
c     -------------------------------------------------------
c     Calculate membrane immersion parameters of a 3D structure
c     given its subunits, borders of each subunit, and
c     its secondary structures (SS)
c
      parameter (maxat=150000,maxres=20000,maxlen=50,
     * maxsegm=900,maxsub=60)
c
c**   'maxat' must be the same as 'maxs' in solva.f !
c
      character*80 pdbtempl,pdbout,line
      character*1 namsu(maxat),namsegm(maxsegm),namesu(maxsub)
      character*4 namat(maxat),namres(maxat)
      character*3 itopo
c
      integer numres(maxat),ienvat(maxat),natsegm(maxsegm),
     * ifirst(maxres),ilast(maxres),ienv(maxres),isegm(2,maxsegm)
c
      real xyz(3,maxat),solv(maxat),xyzbuf(3,20),
     *  xyzc(3,maxat),xyzca(3,maxsegm,maxlen),vcenter(3)
c
      data pi/3.14159/
c
      open (21,file=pdbtempl)
c
c     1. Read atoms from PDB file as a continuous list; define their
c     ASA, solvation parameters, and residue IDs; and prepare
c     list of CA atoms of all SS
c
      read (21,'(a)') line
      if(line(1:21).eq.'REMARK Small molecule') then
         model=1
         call read_small (nat,xyz,solv,numres,namsu,
     *     iprint,namat,namres,nres,ifirst,ilast)
      else
         close(21)
         open (21,file=pdbtempl)
         call readpdb (nat,xyz,solv,numres,namsu,
     *     iprint,namat,namres,nres,ifirst,ilast,ihetero)
      end if
      nsuper=10
      if(nat.lt.10) nsuper=nat
      do i=1,nsuper
         do j=1,3
            xyzbuf(j,i)=xyz(j,i)
         end do
      end do
c
c     Read rotamer library
c
c     Side chain reconstruction and initial optimization
c     (set of atoms will be refedined)
c
c     Find optimal initial orientation and hydrophobic thickness
c     in a hydrophobic slab using method of Tusnady
c
      call orient (nat,xyz,solv,phimax,tetamax,vcenter,
     *  dmax,zcenter,emin)
c
      do i=1,nat
         xyz(1,i)=xyz(1,i)-vcenter(1)
         xyz(2,i)=xyz(2,i)-vcenter(2)
         xyz(3,i)=xyz(3,i)-vcenter(3)
      end do
      phimax=phimax*pi/180.
      tetamax=tetamax*pi/180.
      call tilting (phimax,tetamax,nat,xyz,xyzc)
      do i=1,nat
         xyz(1,i)=xyzc(1,i)
         xyz(2,i)=xyzc(2,i)
         xyz(3,i)=xyzc(3,i) - zcenter
      end do
      write (*,'(''Initial optimiz. in hydrocarbon slab with '',
     *   ''approximate energy functions:'')') 
      write (*,'(''emin='',f8.1,'' slab thickn.='',f8.1)') emin,dmax
      dmin=dmax
c
c     Find initial TM secondary structures if any
c
      nsegm=0
      call find_segm (nat,xyz,xyzca,numres,namat,namsu,
     *  namres,nsegm,isegm,namsegm,natsegm,dmax,itypess)
c
      model=1
      if(nsegm.ge.1) model=2
      if(dmax.gt.14..and.nsegm.eq.0) model=3
c
      dmatch=30.
      dini=15.
      if(model.eq.2.and.itypess.eq.2) then
         dmatch=23.6
         dini=11.8
      end if
c
      ica=0
      write (*,'(i5,'' possible transmembrane secondary structure''
     *  '' segments'')') 
     *  nsegm
      if(nsegm.ge.1) then
         do i=1,nsegm
            write (*,'(i4,1x,a1,2i5,i4)') i,namsegm(i),
     *        (isegm(j,i),j=1,2),natsegm(i)
c           do j=1,natsegm(i)
c              ica=ica+1
c              write (*,'(''CA '',i5,3f8.3)') ica,
c    *           (xyzca(k,i,j),k=1,3)
c           end do
         end do
      end if
c     go to 999
c
c     2. Move center of mass of CA atoms to the origin of coordinates
c
      if(model.eq.1) write (*,'(''Peripheral protein'')')
      if(model.eq.2) write (*,'(''Possible transmembrane protein,'',
     *  '' checking... '')')
      if(model.eq.3)write (*,'(''Possible transmembrane without'',
     *  '' TM secondary structure, checking...'')')
      if (model.eq.1.or.model.eq.3) then
         call bundle_axis0 (nat,namat,xyz)
      else 
         call bundle_axis (nat,xyz,xyzc,nsegm,isegm,natsegm,xyzca,
     *     iprint,numres,namsu,namat,namres)
      end if
c
c     3. Determine "water-facing" atoms and residues.  
c     ASA of water-facing atoms are nullified. Currently,
c     this is not done. 
c
      if(model.eq.1.or.model.eq.3) then
         call watface0 (nat,ienvat,nres,ienv)
      else 
         call watface (nat,xyz,numres,namsu,nres,ifirst,ilast,
     *     ienv,iprint,pdbtempl,nsegm,natsegm,xyzca,
     *     model,ienvat,namat,namres,nlipid)
         write (*,'(''N interact. lipids:'',i4)') nlipid 
      end if
c
      ifunc=2
c
c     Iterative spatial and side-chain optimization
c
c     4. Determine optimal position and hydrophobic thickness
c     immersion depth of the molecule
c
      call optim (model,dmin,dmatch,dini,nsegm,nat,xyz,xyzc,
     *  solv,numres,namsu,shift,tilt,ener,iprint,namat,namres,
     *  pdbtempl,ifunc,d12,t12,nlipid)
c**
c     Define list of names of subunits:
c
      nsu=1
      namesu(1)=namsu(1)
      do i=2,nat
         if(namsu(i).ne.namsu(i-1)) then
            nsu=nsu+1
            namesu(nsu)=namsu(i)
         end if
      end do
c     do i=1,nsu
c        write (*,'(a1)') namesu(i)
c     end do
c
c     Output of tilt angles for individual subunits calculated
c     based on their inertia axes and lists of residues
c     penetratin to the hydrocarbon core
c
      do i=1,nsu
         call depth_angle2 (nat,xyz,tiltot,dmin,dpenetr,
     *     namat,numres,namsu,namesu(i),pdbtempl)
      end do
c**
c     Determine transmembrane segments and re-define their
c     tilt angles and tilt angle of the protein
c
      call depth_angle (nat,xyz,tiltot,dmin,dpenetr)
      if(dpenetr.le.18.) then
         model=1
         dmm=dpenetr+0.4
      else 
         if(nsegm.ge.1) call deftm (nres,ifirst,ilast,nsegm,
     *     namsegm,isegm,nat,xyz,dmin,numres,namsu,iprint,namat,
     *     xyzca,ienvat,solv,pdbtempl,tiltot)
         dmm=dmin+0.8
      end if
c
      write (*,'(a10,1x,''emin='',f8.1,'' thickn='',f5.1,
     *''+-'',f4.1,'' tilt='',f7.0''+-'',f6.0)') 
     *pdbtempl(1:10),ener,dmm,d12,tiltot,t12
      itilt=int(tiltot)
      it12=int(t12)
      write (24,'(a14,'';'',2(f4.1,'';''),2(i4,'';''),
     *  f6.1,'';'')') pdbtempl(1:14),dmm,d12,itilt,it12,ener
c
c     make rotation that provides optimal XY projection for picture
c     of the protein 
c
      call project (nat,xyz)
c
      close (21)
c
c     6. Determine transformation matrix and transform all atoms
c     from the original PDB file including HETERO atoms
c
      open(21,file=pdbtempl)
      open(22,file=pdbout)
c
      call write_pdb (nsuper,xyz,xyzbuf,dmin,itopo,model)
c
      close (21)
      close (22)
c
      return
      end
c
      subroutine optim (model,dmin,dmatch,dini,nsegm,nat,xyz,xyzc,
     *  solv,numres,namsu,shift,tilt,ener,iprint,namat,namres,
     *  pdbtempl,ifunc,d12,t12,nlipid)
c     ------------------------------------------------------------
c     Minimize transfer energy. Free variables for optimization:
c     membrane thickness for TM version; shift of center
c     of mass along Z axis; and rotation angles arond OX and OY axes.
c
      character*80 pdbtempl
      character*1 namsu(1)
      character*4 namat(1),namres(1)
      integer numres(1)
      real xyz(3,1),solv(1),xyzc(3,1),dtab(500),etab(500)
c
      data ecut_tab/1.0/,emism/0.0288/
c     data ecut_tab/1.0/,emism/0.0/
c
      if(model.eq.1) then
         dmembr=30.
         call profile (dmembr)
         call locate (nat,xyz,xyzc,solv,numres,namsu,dmembr,
     *     shift,tilt,ener,iprint,namat,namres,
     *     pdbtempl,ifunc,s12,t12,1)
         dmin=dmembr
         d12=s12
      else
         emin=99.
         dmin=99.
         m=0
         thick1=dini-2.5
         if(model.eq.3.or.(model.eq.2.and.nsegm.le.2)) thick1=10.
         thick2=dini+5.0
         do thick=thick1,thick2,0.20
            dmembr=thick*2.
            call profile (dmembr)
            call locate (nat,xyz,xyzc,solv,numres,namsu,dmembr,
     *        shift,tilt,ener,iprint,namat,namres,
     *        pdbtempl,ifunc,s12,t12,0)
            diff=abs(dmembr-dmatch)
            if(diff.gt.0.5) then
               dmism=diff-0.5
            else
               dmism=0.
            end if
c           write (*,'(2f9.3)') emism,dmism
            ener=ener+float(nlipid)*emism*dmism
            m=m+1
            etab(m)=ener
            dtab(m)=dmembr
c           write (*,'(2f9.3,f8.1)') dtab(m),etab(m),tilt
            if(ener.lt.emin) then
               emin=ener
               dmin=dmembr
            end if
         end do
         ntab=m
         dmin0=100.
         dmax0=0.
         do i=1,ntab
c           write (*,'(2f9.3)') dtab(i),etab(i)
            if(etab(i).le.emin+ecut_tab) then
               if(dtab(i).lt.dmin0) dmin0=dtab(i)
               if(dtab(i).gt.dmax0) dmax0=dtab(i)
            end if
         end do
         d12=0.5*(dmax0-dmin0)
         call profile (dmin)
         call locate (nat,xyz,xyzc,solv,numres,namsu,dmin,
     *     shift,tilt,ener,iprint,namat,namres,
     *     pdbtempl,ifunc,s12,t12,1)
      end if
      return
      end
c
      subroutine write_pdb (nsuper,xyz,xyzbuf,dmin,itopo,model)
c     --------------------------------------------------------
c
      real xyz(3,1),xyzbuf(3,1),vshift1(3),vshift2(3),rot(3,3)
c
      character*3 resname(20)/'LEU','ILE','VAL','PHE','TRP',
     * 'MET','PRO','CYS','TYR','THR','SER','LYS','GLU',
     *'ASP','HIS','GLN','ASN','ARG','ALA','GLY'/,itopo
      character*86 line
c
c     Determine transformation matrix, transform all atoms
c     of the original PDB file including HETERO atoms,
c     and rewite PDB file keeping all other records from
c     original PDB file
c
      npoint=40
      maxpnt=2*npoint+1
c
c     determine transformation matrix and translation vector
c     for moving 'xyzbuf' to the coordinate system of 'xyz' set
c
      call rmsd (nsuper,xyz,xyzbuf,vshift1,vshift2,rot,rms)
c
c     modify coordinates of the original PDB file accordingly:
c
      dmem=dmin/2.
      dmm=dmem+0.4
      write (22,'(''REMARK      1/2 of bilayer thickness:'',f7.1)') dmm
      dxymax=0.
      zav=0.
      topo=1.
      mcount=1
      do i=1,150000
         read (21,'(a)',end=20) line
         if(line(1:6).eq.'ENDMDL') go to 20
         if(line(14:14).eq.'H'.or.line(14:14).eq.'D'.or.
     *     line(13:13).eq.'H'.or.line(14:14).eq.'Q') go to 16
         if(line(18:20).ne.'DUM'.and.line(9:14).ne.
     *      '1/2 of') then
         if(line(1:4).eq.'ATOM'.or.line(1:6).eq.'HETATM') then
            read (line,'(30x,3f8.3)') xc,yc,zc
            xc=xc-vshift2(1)
            yc=yc-vshift2(2)
            zc=zc-vshift2(3)
            a=xc*rot(1,1)+yc*rot(1,2)+zc*rot(1,3)
            b=xc*rot(2,1)+yc*rot(2,2)+zc*rot(2,3)
            c=xc*rot(3,1)+yc*rot(3,2)+zc*rot(3,3)
            a = a + vshift1(1)
            b = b + vshift1(2)
            c = c + vshift1(3)
            if(mcount.eq.1) then
               mcount=2
               if(itopo.eq.'in '.and.c.lt.0.) topo=1.
               if(itopo.eq.'in '.and.c.ge.0.) topo=-1.
               if(itopo.eq.'out'.and.c.lt.0.) topo=-1.
               if(itopo.eq.'out'.and.c.ge.0.) topo=1.
c              write (*,'(''topo='',f4.0)') topo
            end if
            c=c*topo
            b=b*topo
            dxy=sqrt(a*a+b*b)
            if(abs(c).le.25..and.dxy.gt.dxymax) dxymax=dxy
c
            if(line(1:4).eq.'ATOM') then
               do j=1,20
                  if(line(18:20).eq.resname(j)) then
                     write (22,'(a30,3f8.3,a26)') line(1:30),
     *                 a,b,c,line(55:86)
                     go to 15
                  end if
               end do
               write (22,'(''HETATM'',a24,3f8.3,a26)') line(7:30),
     *           a,b,c,line(55:86)
   15          continue
            else
               write (22,'(a30,3f8.3,a26)') line(1:30),
     *           a,b,c,line(55:86)
            end if
            zav=zav+c
         else
           if(line(1:3).ne.'END'.and.line(1:3).ne.'TER'.and.
     *       line(1:6).ne.'ANISOU'.and.line(1:4).ne.'CONE')
     *       write (22,'(a)') line
         end if
         end if
   16    continue
      end do
   20 continue
      m=0
      fn=2.*float(npoint)
      do i=1,maxpnt
         do j=1,maxpnt
            m=m+1
            a=float(2*i)-fn
            b=float(2*j)-fn
            dxy=sqrt(a*a+b*b)
            if(dxy.le.dxymax) then
c              if(.not.(model.eq.1.and.itopo.eq.'out')) then
                  c=-dmem-0.4
                  write (22,'(''HETATM'',i5,2x,''N   '',''DUM '',
     *              1x,i4,4x,3f8.3)') m,m,a,b,c
c              end if
c              if(.not.(model.eq.1.and.itopo.eq.'in ')) then
                  c=dmem+0.4
                  write (22,'(''HETATM'',i5,2x,''O   '',''DUM '',
     *              1x,i4,4x,3f8.3)') m,m,a,b,c
c              end if
            end if
         end do
      end do
      do i=1,150000
         read (21,'(a)',end=30,err=30) line
         if(i.eq.1) write (22,'(''ENDMDL'')') 
         write (22,'(a)') line
      end do
   30 continue
      return
      end
