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
!
!     Format Statements
!
 1000 Format(1x,'Enter Program excitationAnalysis.')
 1010 Format(3x,'Initial State Matrix File: ',A,/,  &
        3x,'Final   State Matrix File: ',A,/)
 2000 Format(/,1x,'***',A,'***',/,  &
        3x,'Number of electrons                            = ',I6,/,  &
        3x,'Number of alpha electrons                      = ',I6,/,  &
        3x,'Number of beta  electrons                      = ',I6,/,  &
        3x,'Number of atomic orbital basis functions       = ',I6,/,  &
        3x,'Number of linearly independent basis functions = ',I6,/)
 3000 Format(/,1x,'Gill Excitation Number: ',F6.2,/,  &
        5x,'Alpha Contribution: ',F6.2,/,  &
        5x,'Beta  Contribution: ',F6.2)
 8999 Format(/,1x,'END OF PROGRAM EXCITATIONANALYSIS.')
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
      nElectronsFinal = GMatrixFileFinal%getVal('nelectrons')
      nElectronsAlphaFinal = GMatrixFileFinal%getVal('nalpha')
      nElectronsBetaFinal = GMatrixFileFinal%getVal('nbeta')
      nBasisFinal = GMatrixFileFinal%getVal('nbasis')
      nBasisUseFinal = GMatrixFileFinal%getVal('nbasisUse')
      ALLOCATE(CAlphaFinal(nBasisFinal,nBasisUseFinal),  &
        CBetaFinal(nBasisFinal,nBasisUseFinal))
      call GMatrixFileFinal%getArray('ALPHA MO COEFFICIENTS',mqcVarOut=tmpMQCvar)
      CAlphaFinal = tmpMQCvar
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
      write(IOut,3000) nExcitationGill,nExcitationAlphaGill,nExcitationBetaGill
!
  999 Continue
      write(iOut,8999)
      end program excitationAnalysis
