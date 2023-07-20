
def get_import_mtsat_cmd(wildcards,input,output):
    t1w_nii = sorted(glob( input.nii_folder + f'/*_MT_T1.nii'))[-1] #get last one only
    mtpdw_nii = sorted(glob( input.nii_folder + f'/*_MT_ISO250_Marm.nii'))[-1] #get last one only
    t1w_prefix=t1w_nii[:-4]
    mtpdw_prefix=mtpdw_nii[:-4]

    cmds=[]

    cmds.append(f'fslsplit {mtpdw_nii} split_mtpd -t')
    cmds.append(f'mtw=`ls split_mtpd* | head -n 1`')
    cmds.append(f'pdw=`ls split_mtpd* | tail -n 1`')
    cmds.append(f'c3d $mtw -orient RAI -o {output.mtw_nii}')
    cmds.append(f'c3d $pdw -orient RAI -o {output.pdw_nii}')

    #MT T1 has reordered axes for some reason (arrghh..)
    cmds.append(f'fslswapdim {t1w_nii} x -z -y {output.t1w_nii}')
    cmds.append(f'fslcpgeom {output.mtw_nii} {output.t1w_nii}')


    cmds.append(f'cp {t1w_prefix}.json {output.t1w_json}')
    cmds.append(f'cp {mtpdw_prefix}.json {output.pdw_json}')
    cmds.append(f'cp {mtpdw_prefix}.json {output.mtw_json}')
    return ' && '.join(cmds)






rule import_mtsat_files:
    input:
        nii_folder='niftis/sub-{subject}',
    params:
        cmd=get_import_mtsat_cmd
    output:
        t1w_nii=bids(root='bids',
                subject='{subject}',
                acq='T1w',
                datatype='anat',
                suffix='MTsat.nii.gz'),
        pdw_nii=bids(root='bids',
                subject='{subject}',
                acq='PDw',
                datatype='anat',
                suffix='MTsat.nii.gz'),
        mtw_nii=bids(root='bids',
                subject='{subject}',
                acq='MTw',
                datatype='anat',
                suffix='MTsat.nii.gz'),
        t1w_json=bids(root='bids',
                subject='{subject}',
                acq='T1w',
                datatype='anat',
                suffix='MTsat.json'),
        pdw_json=bids(root='bids',
                subject='{subject}',
                acq='PDw',
                datatype='anat',
                suffix='MTsat.json'),
        mtw_json=bids(root='bids',
                subject='{subject}',
                acq='MTw',
                datatype='anat',
                suffix='MTsat.json'),
    shadow: 'minimal'
    shell:
        "{params.cmd}"
 

       

