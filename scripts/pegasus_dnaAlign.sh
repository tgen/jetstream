#!/usr/bin/env bash
#####################################################################
# Copyright (c) 2011 by The Translational Genomics Research
# Institute. All rights reserved. This License is limited to, and you may
# use the Software solely for, your own internal and non-commercial use
# for academic and research purposes. Without limiting the foregoing, you
# may not use the Software as part of, or in any way in connection with
# the production, marketing, sale or support of any commercial product or
# service or for any governmental purposes. For commercial or governmental
# use, please contact dcraig@tgen.org. By installing this Software you are
# agreeing to the terms of the LICENSE file distributed with this
# software.
#####################################################################

thisStep="pegasus_nextJob_dnaAlign.txt"
nxtStep1="pegasus_nextJob_indelRealign.txt"
nxtStep2="pegasus_nextJob_recalibrate.txt"

constants=${JETSTREAM_HOME}/centralPipe/constants/constants.txt
constantsDir=${JETSTREAM_HOME}/centralPipe/constants/
myName=`basename $0 | cut -d_ -f2`

time=`date +%d-%m-%Y-%H-%M`
echo "Starting $0 at $time"
if [ "$1" == "" ] ; then
    echo "### Please provide runfolder as the only parameter"
    echo "### Exiting!!!"
    exit
fi
runDir=$1
projName=`basename $runDir | awk -F'_ps20' '{print $1}'`
configFile=$runDir/$projName.config
if [ ! -e $configFile ] ; then
    echo "### Config file not found at $configFile!!!"
    echo "### Exiting!!!"
    exit
else
    echo "### Config file found."
fi
recipe=`cat $configFile | grep "^RECIPE=" | cut -d= -f2 | head -1 | tr -d [:space:]`
debit=`cat $configFile | grep "^DEBIT=" | cut -d= -f2 | head -1 | tr -d [:space:]`

# This was set by recipes, but it is 8 in every single recipe. Bumped to 12
# so that there are cores for samtools to use (bwa uses 8). 
nCores=12


dnaAligner=`grep "@@"$recipe"@@" $constants | grep @@DNAALIGNER= | cut -d= -f2`
ref=`grep "@@"$recipe"@@" $constants | grep @@REF= | cut -d= -f2`
samtoolsPath=`grep "@@"$recipe"@@" $constants | grep @@SAMTOOLSPATH= | cut -d= -f2`
bwaPath=`grep "@@"$recipe"@@" $constants | grep @@BWAPATH= | cut -d= -f2`
faiFile="$ref.fai"

echo "### projName: $projName"
echo "### confFile: $configFile"
echo "### dnaAlign: $dnaAligner"
d=`echo $runDir | cut -c 2-`

skipLines=1
skipGroup=0
qsubFails=0
###
for configLine in `cat $configFile`
do
    if [ "$configLine" == "=START" ] ; then
        skipLines=0
        continue
    fi
    if [ $skipLines -eq 0 ] ; then
        if [[ $configLine == SAMPLE* || $configLine == =END* ]] ; then
            #echo "config line is $configLine"
            arrayCount=${#mergeArray[@]}
            if [ $arrayCount -gt 0 ] ; then
                ((sampleCount++))
                echo "### What I have: Kit: $kitName, sample: $samName, assay: $assayID, libraID: $libraID"
                if [ "$assayID" == "RNA" ] ; then
                    echo "### Assay ID is $assayID. Will skip."
                    skipGroup=1
                fi
                lastItemIndex=`echo ${#mergeArray[@]}-1 | bc`
                if [ $skipGroup -ne 1 ] ; then
                    for (( i=0; i<=$lastItemIndex; i++ ))
                    do
                        #actually this loop is all unncessary... gotta merge some other way
                        #echo "array with index $i is::: ${mergeArray[$i]}"
                        bamPre="$runDir/$kitName/$samName/$samName.preMerge."`printf "%03d" "$i"`
                        bamName="$runDir/$kitName/$samName/$samName.preMerge."`printf "%03d" "$i"`".bam"
                        targetName="$runDir/$kitName/$samName/splitFastqs/$samName"`printf "_%03d" "$i"`"_R1.fastq.gz"
                        echo "### BAMPRE: $bamPre"
                        #rdGrpID=`echo ${mergeArray[$i]} | cut -d= -f2 | cut -d, -f1`
                        #thisFq=`echo ${mergeArray[$i]} | cut -d= -f2 | cut -d, -f2`
                        thisFq="$runDir/$kitName/$samName/$samName"`printf "_%03d" "$i"`"_R1.fastq.gz"
                        thisR2=`echo $thisFq | sed 's/\(.*\)_R1/\1_R2/'`
                        RG_CN=TGen
                        RG_PL=ILLUMINA
                        RG_ID_fromConfig=`echo ${mergeArray[$i]} | cut -d= -f2 | cut -d, -f1`
                        FCID=`zcat $thisFq | head -n1 | cut -d: -f3`
                        LANE=`zcat $thisFq | head -n1 | cut -d: -f4`
                        INDEX=`zcat $thisFq | head -n1 | cut -d: -f10`
                        RG_LB=`echo $RG_ID_fromConfig | cut -d_ -f3`
                        RG_PU="${FCID}_${LANE}"
                        RG_ID="${FCID}_${LANE}_${RG_LB}"

                        if [[ ! -e $thisFq || ! -e $thisR2 ]] ; then
                            echo "### Fastq R1 or R2 itself doesnt exist. Weird: $thisFq"
                            ((qsubFails++))
                            continue
                        fi
                        if [ -e $thisFq.needsToBeSplit ] ; then
                            echo "### This fastq is too large"
                            #marking this file so indel realign or recalibration does not attempt to run it in whole
                            echo "### Marking this file to process split"
                            touch $bamName.processSplit
                            continue
                        fi
                        if [[ -e $bamName.dnaAlignPass || -e $bamName.dnaAlignInQueue || -e $bamName.dnaAlignFail ]] ; then
                            echo "### This bam file already passed, failed, or inQueue."
                            continue
                        fi
                        d=`echo $runDir | cut -c 2-`
                        echo "### Aligner for recipe $recipe is $dnaAligner"
                        case $dnaAligner in
                        bwamem) echo "### Submitting to bwaMem to create $bamName"


                            if [ -z "$INDEX" ] ; then
                                echo "This sample doesn't have an index indicated in the FASTQ"
                                rgTag="@RG\tID:${RG_ID}\tSM:$samName\tPL:${RG_PL}\tCN:${RG_CN}\tPU:${RG_PU}\tLB:${RG_LB}"
                            else
                                rgTag="@RG\tID:${RG_ID}\tSM:$samName\tPL:${RG_PL}\tCN:${RG_CN}\tPU:${RG_PU}\tLB:${RG_LB}\tKS:${INDEX}"
                            fi

                            echo "$rgTag"
                            sbatch --account ${debit} --output $runDir/oeFiles/%x-slurm-%j.out -n 1 -N 1 --cpus-per-task $nCores  --mem-per-cpu 3000 --export ALL,D=$d,RGTAG="$rgTag",FASTQ1=$thisFq,FASTQ2=$thisR2,REF=$ref,BWAPATH=$bwaPath,SAMTOOLSPATH=$samtoolsPath,FAI=$faiFile,BAMPRE=$bamPre,RUNDIR=$runDir,NXT1=$nxtStep1,NXT2=$nxtStep2,D=$d ${JETSTREAM_HOME}/pegasusPipe/jobScripts/pegasus_bwaMem.sh
                            if [ $? -eq 0 ] ; then
                                touch $bamName.dnaAlignInQueue
                            else
                                ((qsubFails++))
                            fi
                            sleep 1
                            ;;
                        bwa) echo "bwa case"
                            ;;
                        *) echo "I should not be here"
                            ;;
                        esac
                        #thisFq=${cpFqPassFile/.cpFqPass/}
                    done
                else #else of skipGroup
                    echo "### Skipping group."
                fi #end of skipGroup
            fi
            kitName=`echo $configLine | cut -d= -f2 | cut -d, -f1`
            samName=`echo $configLine | cut -d= -f2 | cut -d, -f2`
            assayID=`echo $configLine | cut -d= -f2 | cut -d, -f3`
            libraID=`echo $configLine | cut -d= -f2 | cut -d, -f4`
            #insSize=`echo $configLine | cut -d= -f2 | cut -d, -f6`
            unset mergeArray
            count=0
            skipGroup=0
            continue
        else #doesnt start with =, add to mergeArray
            mergeArray[$count]=$configLine
            ((count++))
        fi
    else
        continue
    fi
done

if [ $qsubFails -eq 0 ] ; then
#all jobs submitted succesffully, remove this dir from messages
    echo "### I should remove $thisStep from $runDir."
    rm -f $runDir/$thisStep
else
#qsub failed at some point, this runDir must stay in messages
    echo "### Failure in qsub. Not touching $thisStep"
fi

time=`date +%d-%m-%Y-%H-%M`
echo "Ending $0 at $time"
