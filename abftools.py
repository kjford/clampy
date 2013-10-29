#     Abftools
#     A set of programs for operating on axon binary files
#     Copyright (C) 2013 Kevin Ford
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>

#     Contributions
#     loadabf is derived from Matlab source code (abfload.m)
#     Original version by Harald Hentschke (harald.hentschke@uni-tuebingen.de)
#     Extended to abf version 2.0 by Forrest Collman (fcollman@Princeton.edu)

#     Version 0.0
#     works for abf v1.83 and abf v2 episodic and gapfree

import numpy as np
import struct
import os.path

def loadabf(fileName):
    # opens abf files
    # Input: fileName as string
    # Output: data as dict with keys data and si sampling interval in us

    
    verbose=0
    # some constants
    BLOCKSIZE=512

    # check existence of file
    if not os.path.isfile(fileName):
        print('File %s does not exist' % fileName)
        return None
    
    # header and section information
    # format as dict {key:,[bit position,dtype,size in bytes, num to read]
    # for abf ver2
    headPar={
    'fFileSignature':[0,'c',1,4],
    'fFileVersionNumber':[4,'b',1,4],
    'uFileInfoSize':[8,'I',4,1],
    'lActualEpisodes':[12,'I',4,1],
    'uFileStartDate':[16,'I',4,1],
    'lFileStartTime':[20,'I',4,1],
    'uStopwatchTime':[24,'I',4,1],
    'nFileType':[28,'h',2,1],
    'nDataFormat':[30,'h',2,1],
    'nSimultaneousScan':[32,'h',2,1],
    'nCRCEnable':[34,'h',2,1],
    'uFileCRC':[36,'I',4,1],
    'FileGUID':[40,'I',4,1],
    'uCreatorVersion':[56,'I',4,1],
    'uCreatorNameIndex':[60,'I',4,1],
    'uModifierVersion':[64,'I',4,1],
    'uModifierNameIndex':[68,'I',4,1],
    'uProtocolPathIndex':[72,'I',4,1]
    }
    
    Sections=['ProtocolSection',
    'ADCSection',
    'DACSection',
    'EpochSection',
    'ADCPerDACSection',
    'EpochPerDACSection',
    'UserListSection',
    'StatsRegionSection',
    'MathSection',
    'StringsSection',
    'DataSection',
    'TagSection',
    'ScopeSection',
    'DeltaSection',
    'VoiceTagSection',
    'SyncraySection',
    'AnnotationSection',
    'StatsSection'
    ]
    
    # dict format key:dtype,n values    
    ProtocolInfo=[[
    'nOperationMode','h',2,1],
    ['fADCSequenceInterval','f',4,1],
    ['bEnableFileCompression','x',1,1],
    ['sUnused1','c',1,3],
    ['uFileCompressionRatio','I',4,1],
    ['fSynchTimeUnit','f',4,1],
    ['fSecondsPerRun','f',4,1],
    ['lNumSamplesPerEpisode','i',4,1],
    ['lPreTriggerSamples','i',4,1],
    ['lEpisodesPerRun','i',4,1],
    ['lRunsPerTrial','i',4,1],
    ['lNumberOfTrials','i',4,1],
    ['nAveragingMode','h',2,1],
    ['nUndoRunCount','h',2,1],
    ['nFirstEpisodeInRun','h',2,1],
    ['fTriggerThreshold','f',4,1],
    ['nTriggerSource','h',2,1],
    ['nTriggerAction','h',2,1],
    ['nTriggerPolarity','h',2,1],
    ['fScopeOutputInterval','f',4,1],
    ['fEpisodeStartToStart','f',4,1],
    ['fRunStartToStart','f',4,1],
    ['lAverageCount','i',4,1],
    ['fTrialStartToStart','f',4,1],
    ['nAutoTriggerStrategy','h',2,1],
    ['fFirstRunDelayS','f',4,1],
    ['nChannelStatsStrategy','h',2,1],
    ['lSamplesPerTrace','i',4,1],
    ['lStartDisplayNum','i',4,1],
    ['lFinishDisplayNum','i',4,1],
    ['nShowPNRawData','h',2,1],
    ['fStatisticsPeriod','f',4,1],
    ['lStatisticsMeasurements','i',4,1],
    ['nStatisticsSaveStrategy','h',2,1],
    ['fADCRange','f',4,1],
    ['fDACRange','f',4,1],
    ['lADCResolution','i',4,1],
    ['lDACResolution','i',4,1],
    ['nExperimentType','h',2,1],
    ['nManualInfoStrategy','h',2,1],
    ['nCommentsEnable','h',2,1],
    ['lFileCommentIndex','i',4,1],
    ['nAutoAnalyseEnable','h',2,1],
    ['nSignalType','h',2,1],
    ['nDigitalEnable','h',2,1],
    ['nActiveDACChannel','h',2,1],
    ['nDigitalHolding','h',2,1],
    ['nDigitalInterEpisode','h',2,1],
    ['nDigitalDACChannel','h',2,1],
    ['nDigitalTrainActiveLogic','h',2,1],
    ['nStatsEnable','h',2,1],
    ['nStatisticsClearStrategy','h',2,1],
    ['nLevelHysteresis','h',2,1],
    ['lTimeHysteresis','i',4,1],
    ['nAllowExternalTags','h',2,1],
    ['nAverageAlgorithm','h',2,1],
    ['fAverageWeighting','f',4,1],
    ['nUndoPromptStrategy','h',2,1],
    ['nTrialTriggerSource','h',2,1],
    ['nStatisticsDisplayStrategy','h',2,1],
    ['nExternalTagType','h',2,1],
    ['nScopeTriggerOut','h',2,1],
    ['nLTPType','h',2,1],
    ['nAlternateDACOutputState','h',2,1],
    ['nAlternateDigitalOutputState','h',2,1],
    ['fCellID','f',4,3],
    ['nDigitizerADCs','h',2,1],
    ['nDigitizerDACs','h',2,1],
    ['nDigitizerTotalDigitalOuts','h',2,1],
    ['nDigitizerSynchDigitalOuts','h',2,1],
    ['nDigitizerType','h',2,1]
    ]
    
    EpochInfoPerDAC=[
    ['nEpochNum','h',2,1],
    ['nDACNum','h',2,1],
    ['nEpochType','h',2,1],
    ['fEpochInitLevel','f',4,1],
    ['fEpochLevelInc','f',4,1],
    ['lEpochInitDuration','i',4,1],
    ['lEpochDurationInc','i',4,1],
    ['lEpochPulsePeriod','i',4,1],
    ['lEpochPulseWidth','i',4,1]
    ]


    ADCInfo=[
    ['nADCNum','h',2,1],
    ['nTelegraphEnable','h',2,1],
    ['nTelegraphInstrument','h',2,1],
    ['fTelegraphAdditGain','f',4,1],
    ['fTelegraphFilter','f',4,1],
    ['fTelegraphMembraneCap','f',4,1],
    ['nTelegraphMode','h',2,1],
    ['fTelegraphAccessResistance','f',4,1],
    ['nADCPtoLChannelMap','h',2,1],
    ['nADCSamplingSeq','h',2,1],
    ['fADCProgrammableGain','f',4,1],
    ['fADCDisplayAmplification','f',4,1],
    ['fADCDisplayOffset','f',4,1],
    ['fInstrumentScaleFactor','f',4,1],
    ['fInstrumentOffset','f',4,1],
    ['fSignalGain','f',4,1],
    ['fSignalOffset','f',4,1],
    ['fSignalLowpassFilter','f',4,1],
    ['fSignalHighpassFilter','f',4,1],
    ['nLowpassFilterType','c',1,1],
    ['nHighpassFilterType','c',1,1],
    ['fPostProcessLowpassFilter','f',4,1],
    ['nPostProcessLowpassFilterType','c',1,1],
    ['bEnabledDuringPN','x',1,1],
    ['nStatsChannelPolarity','h',2,1],
    ['lADCChannelNameIndex','i',4,1],
    ['lADCUnitsIndex','i',4,1]
    ]
    
    # for abf ver < 2
    # format as dict {key:,[bit position,dtype,size in bytes, number to read]
    headParold={
    'fFileSignature':[0,'c',1,4],
    'fFileVersionNumber':[4,'f',4,1],
    'nOperationMode':[8,'h',2,1],
    'lActualAcqLength':[10,'i',4,1],
    'nNumPointsIgnored':[14,'h',2,1],
    'lActualEpisodes':[16,'i',4,1],
    'lFileStartTime':[24,'i',4,1],
    'lDataSectionPtr':[40,'i',4,1],
    'lSynchArrayPtr':[92,'i',4,1],
    'lSynchArraySize':[96,'i',4,1],
    'nDataFormat':[100,'h',2,1],
    'nADCNumChannels':[120, 'h',2,1],
    'fADCSampleInterval':[122,'f',4,1],
    'fSynchTimeUnit':[130,'f',4,1],
    'lNumSamplesPerEpisode':[138,'i',4,1],
    'lPreTriggerSamples':[142,'i',4,1],
    'lEpisodesPerRun':[146,'i',4,1],
    'fADCRange':[244, 'f',4,1],
    'lADCResolution':[252, 'i',4,1],
    'nFileStartMillisecs':[366, 'h',2,1],
    'nADCPtoLChannelMap':[378, 'h',2,16],
    'nADCSamplingSeq':[410, 'h',2,16],
    'sADCChannelName':[442, 'B',1,160],
    'fADCProgrammableGain':[730, 'f',4,16],
    'fInstrumentScaleFactor':[922, 'f',4,16],
    'fInstrumentOffset':[986, 'f',4,16],
    'fSignalGain':[1050, 'f',4,16],
    'fSignalOffset':[1114, 'f',4,16],
    'nTelegraphEnable':[4512,'h',2,16],
    'fTelegraphAdditGain':[4576,'f',4,16]
    }

    # open file
    fileO = open(fileName,'rb')
    fileO.seek(0)
    
    # determine version
    verName = fileO.read(4) # this should read in as a str in python 2
    if verName[-1]=='2':
        headParuse = headPar
        verNum=2
    else:
        headParuse = headParold
        verNum=1

    if verbose:
        print('Version: %d' % verNum)
        
    # get header info into H
    H = dict()
    for hkey in headParuse:
        bType = headParuse[hkey][1] * headParuse[hkey][3]
        seekPos = headParuse[hkey][0]
        numBytes = headParuse[hkey][2] * headParuse[hkey][3]
        fileO.seek(seekPos)
        infoStr = fileO.read(numBytes)
        infoByte = struct.unpack(bType,infoStr)
        H[hkey]=(list(infoByte))
    
    # get section info if Ver 2
    # add to the header parameters
    # note: in versions <2 this is contained in the header
    if verNum == 2:
        offset = 76
        secInfo=dict()
        for sec in Sections:
            fileO.seek(offset)
            blockInd = list(struct.unpack('I',fileO.read(4)))
            nBytes = list(struct.unpack('I',fileO.read(4)))
            llEntries = list(struct.unpack('q',fileO.read(8)))
            secInfo[sec]= [blockInd,nBytes,llEntries]
            offset += 16
        # get Protocol infor
        offset = secInfo['ProtocolSection'][0][0] * BLOCKSIZE
        for pkey in ProtocolInfo:
            # add to header dict
            numBytes = pkey[2] * pkey[3]
            bType = pkey[1] * pkey[3]
            fileO.seek(offset)
            infoByte = list(struct.unpack(bType,fileO.read(numBytes)))
            H[pkey[0]] = infoByte
            offset += numBytes
        # might do epoch section later
        nChan = secInfo['ADCSection'][2]
        for chanInd in range(nChan[0]):
            offset = secInfo['ADCSection'][0][0] * BLOCKSIZE + chanInd * secInfo['ADCSection'][2][0]
            for akey in ADCInfo:
                # add to header dict
                numBytes = akey[2] * akey[3]
                bType = akey[1] * akey[3]
                fileO.seek(offset)
                infoByte = list(struct.unpack(bType,fileO.read(numBytes)))
                if chanInd > 0:
                    H[akey[0]].append(infoByte)
                else:
                    H[akey[0]] = infoByte
                offset += numBytes
        H['nADCNumChannels']=secInfo['ADCSection'][2]
        H['lActualAcqLength']=secInfo['DataSection'][2]
        H['lDataSectionPtr']=secInfo['DataSection'][0]
        H['nNumPointsIgnored']=[0]
        
        H['fADCSampleInterval']=[H['fADCSequenceInterval'][0]/H['nADCNumChannels'][0]]
        # headParuse should have everything we need now
    # end if
    
    # Read data
    # ver 0.0 only gapfree and episodic
    if H['nDataFormat'][0] == 0:
        dataSz=2;  # bytes/point
        precision='h'
    elif H['nDataFormat'][0] == 1:
        dataSz=4;  # bytes/point
        precision='f'
    else:
        print('Unknown data format')
        fileO.close()
        return None
    
    nSweeps=H['lActualEpisodes'][0]
    sweeps = np.arange(nSweeps)
    nCh = H['nADCNumChannels'][0]
    chans = H['nADCSamplingSeq'][:nCh]
    si=H['fADCSampleInterval'][0] * nCh
    addGain = []
    if verbose:
        print('File version: %f8' % H['fFileVersionNumber'][-1])
    if H['fFileVersionNumber'][-1]>=1.65:
        for ch in chans:
            chAddGain = H['nTelegraphEnable'][ch] * H['fTelegraphAdditGain'][ch]
            if chAddGain==0:
                addGain.append(1)
            else:
                addGain.append(chAddGain)    
        
    else:
        for ch in chans:
            addGain.append(1)
    
    headOffset=H['lDataSectionPtr'][0] * BLOCKSIZE + H['nNumPointsIgnored'][0] * dataSz
    if verbose:
        print('#Sweeps: %s, si: %s, #chan: %s, precision: %s' % (nSweeps,si,nCh,precision))
        
    if H['nOperationMode'][0] == 5:        # episodic
        if verbose:
            print('episodic mode')
        # determine how many data points
        dataPts=H['lActualAcqLength'][0]
        dataPtsPerChan=dataPts/H['nADCNumChannels'][0]
        dataPtsPerChanPerSweep=dataPtsPerChan/H['lActualEpisodes'][0]
        dataPtsPerSweep=dataPtsPerChanPerSweep * nCh
        # creating numpy array for data:
        # rows are data, columns are ea sweep, third dimension is each channel
        D=np.zeros((dataPtsPerChanPerSweep,nSweeps,nCh))
        
        selectedSegStartInPts=(sweeps * dataPtsPerSweep) * dataSz + headOffset
        # populate D
        for sw in sweeps:
            fileO.seek(selectedSegStartInPts[sw])
            tmpStr = fileO.read(dataPtsPerSweep * dataSz)
            tmpFmt = np.array(struct.unpack(precision * dataPtsPerSweep,tmpStr))
            tmpbychanint = tmpFmt.reshape(dataPtsPerChanPerSweep,nCh)
            tmpbychan = tmpbychanint.astype('float')    # to do math
            chcount=0
            for ch in chans:
                if H['nDataFormat'][0]==0:
                    tmpbychan[:,chcount] = tmpbychan[:,chcount]/(H['fInstrumentScaleFactor'][ch]* \
                    H['fSignalGain'][ch]*H['fADCProgrammableGain'][ch] *addGain[chcount]) \
                    * H['fADCRange'][0] / H['lADCResolution'][0] \
                    +H['fInstrumentOffset'][ch]-H['fSignalOffset'][ch]
                
                D[:,sw,chcount]=tmpbychan[:,chcount]
                chcount+=1
       
        
    elif H['nOperationMode'][0] == 3:    # gapfree: note: reads beginning to end
        if verbose:
            print('gapfree')    
        dataPtsPerChan=H['lActualAcqLength'][0]/H['nADCNumChannels'][0]
        dataPts=dataPtsPerChan*nCh
        # creating numpy array for data:
        # rows are data, columns are each channel
        D=np.zeros((dataPtsPerChan,nCh))
        fileO.seek(headOffset)
        tmpStr = fileO.read(dataPts * dataSz)
        tmpFmt = np.array(struct.unpack(precision * dataPts,tmpStr))
        tmpbychanint = tmpFmt.reshape(dataPtsPerChan,nCh)
        tmpbychan = tmpbychanint.astype('float')    # to do math
        chcount=0
        for ch in chans:
            if H['nDataFormat'][0]==0:
                tmpbychan[:,chcount] = tmpbychan[:,chcount]/(H['fInstrumentScaleFactor'][ch]* \
                H['fSignalGain'][ch]*H['fADCProgrammableGain'][ch] *addGain[chcount]) \
                * H['fADCRange'][0] / H['lADCResolution'][0] \
               +H['fInstrumentOffset'][ch]-H['fSignalOffset'][ch]
               
            D[:,chcount]=tmpbychan[:,chcount]
            chcount+=1
    else:
        if verobse:
            fileO.close()
            print('Recording mode not supported')
            return None
    fileO.close()
    data = {'Data':D,'si':si}
    return data    
# end loadabf
