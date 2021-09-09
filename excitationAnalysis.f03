INCLUDE 'excitationAnalysis_mod.f03'
      Program excitationAnalysis
!
!     This program uses Gaussian matrix files to carry out various excitation
!     anlysis models.
!
!     -H. P. Hratchian, 2021.
!
!
!     USE Connections
!
      use excitationAnalysis_mod
!
!     Variable Declarations
!
      implicit none
      integer(kind=int64)::nCommands,nElectronsInitial,  &
        nElectronsAlphaInitial,nElectronsBetaInitial,  &
        nElectronsFinal,nElectronsAlphaFinal,nElectronsBetaFinal,  &
        nBasisInitial,nBasisUseInitial,nBasisFinal,nBasisUseFinal
      real(kind=real64)::nExcitationAlphaGill,nExcitationBetaGill,  &
        nExcitationGill
      real(kind=real64),dimension(:,:),allocatable::OverlapAO,  &
        CAlphaInitial,CBetaInitial,CAlphaFinal,CBetaFinal,bMatrixAlpha,  &
        bMatrixBeta
      character(len=512)::matrixFilenameInitial,matrixFilenameFinal
      type(mqc_gaussian_unformatted_matrix_file)::GMatrixFileInitial,  &
        GMatrixFileFinal
      type(MQC_Variable)::tmpMQCvar

!hph+
!      integer(kind=int64)::nCommands,i,j,k1,k2,nAtoms,nAt3
!      integer(kind=int64),dimension(:),allocatable::atomicNumbers
!      real(kind=real64)::Vnn,Escf
!      real(kind=real64),dimension(3)::tmp3Vec
!      real(kind=real64),dimension(:),allocatable::cartesians
!      real(kind=real64),dimension(:,:),allocatable::distanceMatrix
!      character(len=512)::matrixFilename,tmpString
!      type(MQC_Variable)::nEalpha,nEbeta,nEtot,KEnergy,VEnergy,OneElEnergy,  &
!        TwoElEnergy,scfEnergy
!      type(MQC_Variable)::SMatrixAO,TMatrixAO,VMatrixAO,HCoreMatrixAO,  &
!        FMatrixAlpha,FMatrixBeta,PMatrixAlpha,PMatrixBeta,PMatrixTotal,  &
!        ERIs,JMatrixAlpha,KMatrixAlpha
!      type(MQC_R4Tensor)::tmpR4
!hph-

!
!     Format Statements
!
 1000 Format(1x,'Enter Program excitationAnalysis.')
 1010 Format(3x,'Initial State Matrix File: ',A,/,  &
        3x,'Final   State Matrix File: ',A,/)
 2000 Format(1x,'***',A,'***',/,  &
        3x,'Number of electrons                            = ',I6,/,  &
        3x,'Number of alpha electrons                      = ',I6,/,  &
        3x,'Number of beta  electrons                      = ',I6,/,  &
        3x,'Number of atomic orbital basis functions       = ',I6,/,  &
        3x,'Number of linearly independent basis functions = ',I6,/)

      
 1100 Format(1x,'nAtoms=',I4)
 1200 Format(1x,'Atomic Coordinates (Angstrom)')
 1210 Format(3x,I3,2x,A2,5x,F7.4,3x,F7.4,3x,F7.4)
 1300 Format(1x,'Nuclear Repulsion Energy = ',F20.6)
 8999 Format(/,1x,'END OF TEST PROGRAM scfEnergyTerms.')
!
!
      write(IOut,1000)
!
!     Get the command line arguments.

      nCommands = command_argument_count()
      if(nCommands.eq.0)  &
        call mqc_error('No command line arguments provided. The input Gaussian matrix file names are required.')
      call get_command_argument(1,matrixFilenameInitial)
      call get_command_argument(2,matrixFilenameFinal)
      write(IOut,1010) TRIM(matrixFilenameInitial),TRIM(matrixFilenameFinal)
!
!     Open the initial state matrix file and read what we need from that file.
!     Once things are read, close the file.
!
      call GMatrixFileInitial%load(matrixFilenameInitial)
      nElectronsInitial = GMatrixFileInitial%getVal('nelectrons')
      nElectronsAlphaInitial = GMatrixFileInitial%getVal('nalpha')
      nElectronsBetaInitial = GMatrixFileInitial%getVal('nbeta')
      nBasisInitial = GMatrixFileInitial%getVal('nbasis')
      nBasisUseInitial = GMatrixFileInitial%getVal('nbasisUse')
      ALLOCATE(OverlapAO(nBasisInitial,nBasisInitial),  &
        CAlphaInitial(nBasisInitial,nBasisUseInitial),  &
        CBetaInitial(nBasisInitial,nBasisUseInitial))
      call GMatrixFileInitial%getArray('OVERLAP',mqcVarOut=tmpMQCvar)
      OverlapAO = tmpMQCvar
      call GMatrixFileInitial%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=tmpMQCvar)
      CAlphaInitial = tmpMQCvar
      Write(IOut,*)' Initial State isUnrestricted=',GMatrixFileInitial%isUnrestricted()
      if(GMatrixFileInitial%isUnrestricted()) then
        call GMatrixFileInitial%getArray('BETA MO COEFFICIENTS',mqcVarOut=tmpMQCvar)
        CBetaInitial = tmpMQCvar
      else
        CBetaInitial = CAlphaInitial
      endIf
      call GMatrixFileInitial%closeFile()
      write(IOut,2000) 'Initial State',nElectronsInitial,nElectronsAlphaInitial,  &
        nElectronsBetaInitial,nBasisInitial,nBasisUseInitial
!
!     Open the final state matrix file and read what we need from that file.
!     Once things are read, close the file.
!
      call GMatrixFileFinal%load(matrixFilenameFinal)
      write(IOut,*)' Final   ICGU = ',GMatrixFileFinal%ICGU
      nElectronsFinal = GMatrixFileFinal%getVal('nelectrons')
      nElectronsAlphaFinal = GMatrixFileFinal%getVal('nalpha')
      nElectronsBetaFinal = GMatrixFileFinal%getVal('nbeta')
      nBasisFinal = GMatrixFileFinal%getVal('nbasis')
      nBasisUseFinal = GMatrixFileFinal%getVal('nbasisUse')
      ALLOCATE(CAlphaFinal(nBasisFinal,nBasisUseFinal),  &
        CBetaFinal(nBasisFinal,nBasisUseFinal))
      call GMatrixFileFinal%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=tmpMQCvar)
      CAlphaFinal = tmpMQCvar
      Write(IOut,*)' Final   State isUnrestricted=',GMatrixFileFinal%isUnrestricted()
      if(GMatrixFileFinal%isUnrestricted()) then
        call GMatrixFileFinal%getArray('BETA MO COEFFICIENTS',mqcVarOut=tmpMQCvar)
        CBetaFinal = tmpMQCvar
      else
        CBetaFinal = CAlphaFinal
      endIf
      call GMatrixFileFinal%closeFile()
      write(IOut,2000) 'Final State',nElectronsFinal,nElectronsAlphaFinal,  &
        nElectronsBetaFinal,nBasisFinal,nBasisUseFinal
!
!     Allocate space for bMatrix, which gives MO expansion coefficients for
!     the final state's MOs in the initial state's MO basis.
!
      ALLOCATE(bMatrixAlpha(nElectronsBetaInitial,nElectronsAlphaFinal),  &
        bMatrixBeta(nElectronsBetaInitial,nElectronsBetaFinal))
      Write(IOut,*)' Hrant - bMatrixAlpha:  ',SHAPE(bMatrixAlpha)
      Write(IOut,*)' Hrant - CAlphaInitial: ',SHAPE(CAlphaInitial)
      Write(IOut,*)' Hrant - CAlphaFinal:   ',SHAPE(CAlphaFinal)
      bMatrixAlpha = MatMul(Transpose(CAlphaInitial(:,1:nElectronsAlphaInitial)),  &
        MatMul(OverlapAO,CAlphaFinal(:,1:nElectronsAlphaFinal)))
      call MQC_Print(IOut,bMatrixAlpha,Header='Alpha b')
      bMatrixBeta = MatMul(Transpose(CBetaInitial(:,1:nElectronsBetaInitial)),  &
        MatMul(OverlapAO,CBetaFinal(:,1:nElectronsBetaFinal)))
      call MQC_Print(IOut,bMatrixBeta,Header='Beta b')

      nExcitationAlphaGill =  nElectronsAlphaInitial -   &
        dot_product(RESHAPE(bMatrixAlpha,(/nElectronsAlphaInitial*nElectronsAlphaFinal/)),  &
        RESHAPE(bMatrixAlpha,(/nElectronsAlphaInitial*nElectronsAlphaFinal/)))
      nExcitationBetaGill =  nElectronsBetaInitial -   &
        dot_product(RESHAPE(bMatrixBeta,(/nElectronsBetaInitial*nElectronsBetaFinal/)),  &
        RESHAPE(bMatrixBeta,(/nElectronsBetaInitial*nElectronsBetaFinal/)))
      nExcitationGill = nExcitationAlphaGill + nExcitationBetaGill
      Write(IOut,*)' Gill excitation number (alpha) = ',nExcitationAlphaGill
      Write(IOut,*)' Gill excitation number (beta ) = ',nExcitationBetaGill
      Write(IOut,*)' Gill excitation number (total) = ',nExcitationGill
!
  999 Continue
      write(iOut,8999)
      end program excitationAnalysis
