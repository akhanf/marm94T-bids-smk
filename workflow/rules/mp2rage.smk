# grabs inversion images (complex, channel data), and performs mp2rage processing
# TODO:
#  - [x] create UNI-DEN image for each channel
#  - [ ] perform channel combination 
#  - [ ] perform T1 mapping  (qMRlab? pymp2rage?)


#mp2rage is  complex conjugate of inv1 * inv 2
# (a-ib)*(c+id)
#  ac +iad - ibc +bd
#  ac+bd +i(bc-ad)
# Re(inv1)*Re(inv2) + Im(inv1)*Im(inv2)  + i ( Im(inv1) * Re(inv2)  - Im(inv1)*Im(inv2)  )
# numerator = Re(above) = Re(inv1)*Re(inv2) + Im(inv1)*Im(inv2)
# denominator = a^2 + b^2 + c^2 + d^2

def get_avg_mp2rage_cmd(wildcards,input,output):
    niftis = glob( input.nii_folder + f'/*cfmmMP2RAGE_3D_ISO250*{wildcards.part}*.nii')
    cmds = []
    cmds.append('c4d')
    cmds.extend(niftis)
    for i in range(len(niftis)-1):
        cmds.append('-add')
    cmds.append(f'-scale {1.0/len(niftis)}')
    cmds.append(f'-scale {1.0/1000}') #rescale lower (values were exploding)
    cmds.append(f'-o {output}')
    return ' '.join(cmds)


rule avg_mp2rage_complex_channels:
    """ average real or imag scans (no moco here yet, anesthetized anyhow..)"""
    input:
        nii_folder='niftis/sub-{subject}'
    params:
        cmd = get_avg_mp2rage_cmd
    output:
        bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',part='{part}',suffix='inversions.nii.gz')
    shell:
        '{params.cmd}'

rule split_inversions:
    input:
        bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',part='{part}',suffix='inversions.nii.gz')
    params:
        vols_inv1='0:7',
        vols_inv2='8:-1'
    output:
        inv1=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',part='{part}',inv='1',suffix='MP2RAGE.nii.gz'),
        inv2=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',part='{part}',inv='2',suffix='MP2RAGE.nii.gz')
    shell:
        'c4d {input} -slice w {params.vols_inv1} -tile w -o {output.inv1} && '
        'c4d {input} -slice w {params.vols_inv2} -tile w -o {output.inv2} '

rule mp2rage_numerator_denominator:
    input:
        re_inv1=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',part='real',inv='1',suffix='MP2RAGE.nii.gz'),
        re_inv2=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',part='real',inv='2',suffix='MP2RAGE.nii.gz'),
        im_inv1=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',part='imag',inv='1',suffix='MP2RAGE.nii.gz'),
        im_inv2=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',part='imag',inv='2',suffix='MP2RAGE.nii.gz'),
    output:
        numerator=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',suffix='MP2RAGEnumerator.nii.gz'),
        denominator=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',suffix='MP2RAGEdenominator.nii.gz'),
    shell:
        'c4d {input.re_inv1} {input.re_inv2} -multiply '
        ' {input.im_inv1} {input.im_inv2} -multiply '
        ' -add -as NUMERATOR -o {output.numerator} '
        ' {input.re_inv1} -dup -multiply '
        ' {input.im_inv1} -dup -multiply '
        ' {input.re_inv2} -dup -multiply '
        ' {input.im_inv2} -dup -multiply '
        ' -add -add -add -as DENOMINATOR -o {output.denominator} '
   
rule mp2rage_add_bias_term: 
    input:
        numerator=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',suffix='MP2RAGEnumerator.nii.gz'),
        denominator=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',suffix='MP2RAGEdenominator.nii.gz'),
    params:
        num_offset = lambda wildcards: '-shift {offset}'.format(offset=-float(wildcards.beta)),
        den_offset = lambda wildcards: '-shift {offset}'.format(offset=2*float(wildcards.beta)),
    output:
        numerator=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',beta='{beta}',suffix='MP2RAGEnumerator.nii.gz'),
        denominator=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',beta='{beta}',suffix='MP2RAGEdenominator.nii.gz'),
    shell:
        'c4d {input.numerator} {params.num_offset} -o {output.numerator} && '
        'c4d {input.denominator} {params.den_offset} -o {output.denominator} '

rule mp2rage_division:
    input:
        numerator=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',beta='{beta}',suffix='MP2RAGEnumerator.nii.gz'),
        denominator=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',beta='{beta}',suffix='MP2RAGEdenominator.nii.gz'),
    output:
        uni=bids(root='work',subject='{subject}',acquisition='mp2rage',datatype='anat',beta='{beta}',suffix='UNIchannels.nii.gz'),
    shell: 
        #c4d divide was not working properly, so using fsl:
        #'c4d {input.numerator} {input.denominator} -divide -replace inf 1000 -inf -1000 NaN 0  -o {output.uni}' 
        'fslmaths {input.numerator} -div {input.denominator} {output.uni}'
 

