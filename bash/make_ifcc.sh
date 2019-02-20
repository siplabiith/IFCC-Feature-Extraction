#!/bin/bash

# Copyright 2012  Johns Hopkins University (Author: Daniel Povey)
# Apache 2.0
# To be run from .. (one directory up from here)
# see ../run.sh for example

# Begin configuration section.
nj=4
cmd=run.pl
compress=true
# End configuration section.

echo "$0 $@"  # Print the command line for logging

if [ -f path.sh ]; then . ./path.sh; fi
. parse_options.sh || exit 1;

if [ $# != 3 ]; then
   echo "Usage: $0 [options] <data-dir> <log-dir> <path-to-ifccdir>";
   echo "e.g.: $0 data/train exp/make_ifcc/train ifcc"
   echo "options: "
   echo "  --ifcc-config <config-file>                      # config passed to compute-ifcc-feats "
   echo "  --nj <nj>                                        # number of parallel jobs"
   echo "  --cmd (utils/run.pl|utils/queue.pl <queue opts>) # how to run jobs."
   exit 1;
fi

data=$1
logdir=$2
ifccdir=$3


# make $ifccdir an absolute pathname.
ifccdir=`perl -e '($dir,$pwd)= @ARGV; if($dir!~m:^/:) { $dir = "$pwd/$dir"; } print $dir; ' $ifccdir ${PWD}`

# use "name" as part of name of the archive.
name=`basename $data`
mkdir -p $ifccdir || exit 1;
mkdir -p $logdir || exit 1;

if [ -f $data/feats.scp ]; then
  mkdir -p $data/.backup
  echo "$0: moving $data/feats.scp to $data/.backup"
  mv $data/feats.scp $data/.backup
fi

scp=$data/wav.scp

required="$scp $ifcc_config"

for f in $required; do
  if [ ! -f $f ]; then
    echo "make_ifcc.sh: no such file $f"
    exit 1;
  fi
done
utils/validate_data_dir.sh --no-text --no-feats $data || exit 1;

for n in $(seq $nj); do
  # the next command does nothing unless $ifccdir/storage/ exists, see
  # utils/create_data_link.pl for more info.
  utils/create_data_link.pl $ifccdir/raw_ifcc_$name.$n.ark
done


echo "$0: [info]: no segments file exists: assuming wav.scp indexed by utterance."
split_scps=""
for n in $(seq $nj); do
    split_scps="$split_scps $logdir/wav_${name}.$n.scp"
done

utils/split_scp.pl $scp $split_scps || exit 1;


  # add ,p to the input rspecifier so that we can just skip over
  # utterances that have bad wave data.

$cmd JOB=1:$nj $logdir/make_ifcc_${name}.JOB.log \
    compute-ifcc-feats $logdir/wav_${name}.JOB.scp \| \
      copy-feats --compress=$compress ark:- \
      ark,scp:$ifccdir/raw_ifcc_$name.JOB.ark,$ifccdir/raw_ifcc_$name.JOB.scp \
     || exit 1;

if [ -f $logdir/.error.$name ]; then
  echo "Error producing ifcc features for $name:"
  tail $logdir/make_ifcc_${name}.1.log
  exit 1;
fi


# concatenate the .scp files together.
for n in $(seq $nj); do
  cat $ifccdir/raw_ifcc_$name.$n.scp || exit 1;
done > $data/feats.scp

rm $logdir/wav_${name}.*.scp  $logdir/segments.* 2>/dev/null

nf=`cat $data/feats.scp | wc -l`
nu=`cat $data/utt2spk | wc -l`
if [ $nf -ne $nu ]; then
  echo "It seems not all of the feature files were successfully processed ($nf != $nu);"
  echo "consider using utils/fix_data_dir.sh $data"
fi

if [ $nf -lt $[$nu - ($nu/20)] ]; then
  echo "Less than 95% the features were successfully generated.  Probably a serious error."
  exit 1;
fi

echo "Succeeded creating IFCC features for $name"
