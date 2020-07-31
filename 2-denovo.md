
#######################################  
## DRAP assemble _de novo_ transcriptomes
#######################################  
### 1.install DRAP
docker pull sigenae/drap
docker create --name drap --privileged -v /mnt/scratch:/docker/scratch -i -t sigenae/drap:latest /bin/bash
### construct the connection between the directory of fastq files and DRAP in docker
#### -u: user id of your account
### 2. run DRAP

### _A. curacao_
mkdir Acura_cat_total
cd Acura_cat_total
cp ../Acura\*.fastq.gz ./  
cat Acura\*\_1.fastq.gz Acura_total_1.fastq.gz  
cat Acura\*\_2.fastq.gz Acura_total_2.fastq.gz  
export SHARED_DIR=$PWD  
mkdir Drap  
nohup docker run --rm -e LOCAL_USER_ID=`id -u \$USER` \

        -u 1003:1003 -v $SHARED_DIR:$SHARED_DIR \
        -w `pwd` sigenae/drap runDrap -o Drap \
        -1 Acura_total_1.fastq.gz\
        -2 Acura_total_2.fastq.gz\
        -s FR --dbg trinity \
        --dbg-mem 128 -norm-mem 128 --no-trim --no-rate --run --write > runDrap-Acura.process 2>&1 &

### _D. aruanus_  
mkdir Daru_cat_total
cd Daru_cat_total
cp ../Daru\*.fastq.gz ./  
cat Daru\*\_1.fastq.gz Daru_total_1.fastq.gz  
cat Daru\*\_2.fastq.gz Daru_total_2.fastq.gz  
export SHARED_DIR=$PWD  
mkdir Drap  
nohup docker run --rm -e LOCAL_USER_ID=`id -u \$USER` \

        -u 1003:1003 -v $SHARED_DIR:$SHARED_DIR \
        -w `pwd` sigenae/drap runDrap -o Drap \
        -1 Daru_total_1.fastq.gz\
        -2 Daru_total_2.fastq.gz\
        -s FR --dbg trinity \
        --dbg-mem 128 -norm-mem 128 --no-trim --no-rate --run --write > runDrap-Daru.process 2>&1 &

### _A. compressus_
mkdir Ocomp_cat_total
cd Ocomp_cat_total
cp ../Ocomp\*.fastq.gz ./  
cat Ocomp\*\_1.fastq.gz Ocomp_total_1.fastq.gz  
cat Ocomp\*\_2.fastq.gz Ocomp_total_2.fastq.gz  
export SHARED_DIR=$PWD  
mkdir Drap  
nohup docker run --rm -e LOCAL_USER_ID=`id -u \$USER` \

        -u 1003:1003 -v $SHARED_DIR:$SHARED_DIR \
        -w `pwd` sigenae/drap runDrap -o Drap \
        -1 Ocomp_total_1.fastq.gz\
        -2 Ocomp_total_2.fastq.gz\
        -s FR --dbg trinity \
        --dbg-mem 128 -norm-mem 128 --no-trim --no-rate --run --write > runDrap-Ocomp.process 2>&1 &

### _P. adelus_
mkdir Padel_cat_total
cd Padel_cat_total
cp ../Padel\*.fastq.gz ./  
cat Padel\*\_1.fastq.gz Padel_total_1.fastq.gz  
cat Padel\*\_2.fastq.gz Padel_total_2.fastq.gz  
export SHARED_DIR=$PWD  
mkdir Drap  
nohup docker run --rm -e LOCAL_USER_ID=`id -u \$USER` \

        -u 1003:1003 -v $SHARED_DIR:$SHARED_DIR \
        -w `pwd` sigenae/drap runDrap -o Drap \
        -1 Padel_total_1.fastq.gz\
        -2 Padel_total_2.fastq.gz\
        -s FR --dbg trinity \
        --dbg-mem 128 -norm-mem 128 --no-trim --no-rate --run --write > runDrap-Padel.process 2>&1 &

### _P. molluscensis_
mkdir Pmol_cat_total
cd Pmol_cat_total
cp ../Pmol\*.fastq.gz ./  
cat Pmol\*\_1.fastq.gz Pmol_total_1.fastq.gz  
cat Pmol\*\_2.fastq.gz Pmol_total_2.fastq.gz  
export SHARED_DIR=$PWD  
mkdir Drap  
nohup docker run --rm -e LOCAL_USER_ID=`id -u \$USER` \

        -u 1003:1003 -v $SHARED_DIR:$SHARED_DIR \
        -w `pwd` sigenae/drap runDrap -o Drap \
        -1 Pmol_total_1.fastq.gz\
        -2 Pmol_total_2.fastq.gz\
        -s FR --dbg trinity \
        --dbg-mem 128 -norm

### 3. Finish the assemble
all_contigs.second_pass.fa in Drap_trinity/e-rmbt_editing would be selected as the _de novo_ assembled transcriptome
