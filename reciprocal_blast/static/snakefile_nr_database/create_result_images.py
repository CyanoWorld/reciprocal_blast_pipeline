#this script is used by the snakemake pipeline in order to produce result images
import matplotlib.pyplot as plt, mpld3
import pandas as pd
#import numpy as np

rec_prot=pd.read_table(snakemake.input['rec_res'])
fw_res=pd.read_table(snakemake.input['fw_res'],header=None)
fw_res.columns=["qseqid", "sseqid", "evalue", "bitscore", "qgi", "sgi", "sacc", "staxids", "sscinames", "scomnames",
                  "stitle"]

fw_res['qseqid'] = fw_res['qseqid'].map(lambda line: line.split('.')[0])
rec_prot = rec_prot.rename(columns={"forward_genome_id": "sacc"})
rec_prot = rec_prot.rename(columns={"backward_genome_id": "qseqid"})
result_data = rec_prot.merge(fw_res,how='inner', on=['sacc','qseqid'])
result_data = result_data.drop_duplicates('sacc', keep='first')

if len(result_data['qseqid'].unique()) <= 200:

    tax_ids = {}
    for query_id in result_data['qseqid'].unique():
        tax_ids[query_id] = []
        for taxonomic_nodes in result_data[result_data['qseqid'] == query_id]['staxids']:
            for taxid in taxonomic_nodes.split(';'):
                if taxid not in tax_ids[query_id]: tax_ids[query_id].append(taxid)

    query_amount = []
    reciprocal_hits = []
    labels = []
    for key, index in zip(tax_ids.keys(), range(1, len(tax_ids.keys()) + 1)):
        reciprocal_hits.append(len(tax_ids[key]))
        query_amount.append(index)
        labels.append(key)

    figure = plt.figure(figsize=[10, 6])
    bar = plt.bar(query_amount, reciprocal_hits, align='center')
    plt.ylabel("amount of organisms in which orthologs are found", fontsize=14)
    plt.xlabel("protein accession number (gi)", fontsize=14)

    for rect in bar:
        height = rect.get_height()
        # print(height)
        plt.text(rect.get_x() + rect.get_width() / 2.0, height, '%d' % int(height), ha='center', va='bottom', fontsize=12)

    input_query_ids = result_data['qseqid'].unique()

    '''
    try:
        query_seqs = open(snakemake.params['input_queries'])
        description = {}
        for line in query_seqs.readlines():
            if '>' in line:
                identifier = line.split(">")[1].split(" ")[0]
                if "." in identifier:
                    identifier = identifier.split(".")[0]
                description[identifier] = line.split(">")[1]
        query_seqs.close()

        for i, box in enumerate(bar.get_children()):
            tooltip = mpld3.plugins.LineLabelTooltip(box, label=description[input_query_ids[i]])
            mpld3.plugins.connect(figure, tooltip)

    except:
        for i, box in enumerate(bar.get_children()):
            tooltip = mpld3.plugins.LineLabelTooltip(box, label=input_query_ids[i])
            mpld3.plugins.connect(figure, tooltip)
    '''
    plt.grid()
    plt.xticks(query_amount, labels, rotation=90)
    plt.savefig(snakemake.params['org_orth_png'])
    plt.savefig(snakemake.output[0])
    mpld3.save_html(figure,snakemake.output[1])



    query_hits_dict = {}
    for prot_id in input_query_ids:
        query_hits_dict[prot_id] = []
        for bw_prot_id in result_data[result_data['qseqid'] == prot_id]['sacc']:
            query_hits_dict[prot_id].append(bw_prot_id)

    query_amount = list(range(1,len(input_query_ids)+1))
    reciprocal_hits = []
    labels = []
    for prot_id in input_query_ids:
        labels.append(prot_id)
        reciprocal_hits.append(len(query_hits_dict[prot_id]))

    figure = plt.figure(figsize=[10,6])
    bar = plt.bar(query_amount, reciprocal_hits, align='center')
    plt.ylabel("amount of hits in other organisms",fontsize=14)
    plt.xlabel("protein accession number (gi)",fontsize=14)
    '''
    try:
        for i, box in enumerate(bar.get_children()):
            tooltip = mpld3.plugins.LineLabelTooltip(box, label=description[input_query_ids[i]])
            mpld3.plugins.connect(figure, tooltip)

    except:
        for i, box in enumerate(bar.get_children()):
            tooltip = mpld3.plugins.LineLabelTooltip(box, label=input_query_ids[i])
            mpld3.plugins.connect(figure, tooltip)
    '''
    for rect in bar:
        height = rect.get_height()
        #print(height)
        plt.text(rect.get_x() + rect.get_width()/2.0, height, '%d' % int(height), ha='center', va='bottom',fontsize=12)


    plt.grid()
    plt.xticks(query_amount, labels, rotation=90)
    plt.savefig(snakemake.params['amount_hits_png'])
    plt.savefig(snakemake.output[2])
    mpld3.save_html(figure,snakemake.output[3])

elif len(result_data['qseqid'].unique()) > 200:
    input_query_ids = result_data['qseqid'].unique()

    tax_ids = {}
    for query_id in result_data['qseqid'].unique():
        tax_ids[query_id] = []
        for taxonomic_nodes in result_data[result_data['qseqid'] == query_id]['staxids']:
            for taxid in taxonomic_nodes.split(';'):
                if taxid not in tax_ids[query_id]: tax_ids[query_id].append(taxid)

    reciprocal_hits = []
    for key, index in zip(tax_ids.keys(), range(1, len(tax_ids.keys()) + 1)):
        reciprocal_hits.append(len(tax_ids[key]))

    query_hits_dict = {}
    for prot_id in input_query_ids:
        query_hits_dict[prot_id] = []
        for bw_prot_id in result_data[result_data['qseqid'] == prot_id]['sacc']:
            if bw_prot_id not in query_hits_dict[prot_id]:
                query_hits_dict[prot_id].append(bw_prot_id)

    query_amount = list(range(1, len(input_query_ids) + 1))
    query_hits = []
    for prot_id in input_query_ids:
        query_hits.append(len(query_hits_dict[prot_id]))


    hbins = round((len(reciprocal_hits) / 4))
    figure_to_png = plt.figure(figsize=[10, 6])
    hist_to_png = plt.hist(reciprocal_hits, bins=hbins)
    hist_to_png = plt.axvline(sum(reciprocal_hits) / len(reciprocal_hits), color='k', linestyle='dashed', linewidth=2,
                              label='average ' + str(round(sum(reciprocal_hits) / len(reciprocal_hits))))
    hist_to_png = plt.axvline(min(reciprocal_hits), color='r', linestyle='dashed', linewidth=2,
                              label='minimum ' + str(min(reciprocal_hits)))
    hist_to_png = plt.axvline(max(reciprocal_hits), color='green', linestyle='dashed', linewidth=2,
                              label='maximum ' + str(max(reciprocal_hits)))
    #hist_to_png = plt.xticks(np.arange(0, max(reciprocal_hits)))
    hist_to_png = plt.xlabel("amount of organisms that have orthologous sequences")
    hist_to_png = plt.ylabel("amount of query sequences")
    plt.savefig(snakemake.params['org_orth_png'])
    plt.savefig(snakemake.output[0])



    figure_to_html = plt.figure(figsize=[10, 6])
    label = ['average ' + str(
    round(sum(reciprocal_hits) / len(reciprocal_hits))),'minimum ' + str(
    min(reciprocal_hits)), 'maximum ' + str(max(reciprocal_hits))]
    hist_to_html = plt.hist(reciprocal_hits, bins=hbins,label=label)
    hist_to_html = plt.xlabel("amount of organisms that have orthologous sequences")
    hist_to_html = plt.ylabel("amount of query sequences")
    hist_to_html = plt.grid(color='grey', linestyle='dashed')
    hist_to_html = plt.legend()
    mpld3.save_html(figure_to_html, snakemake.output[1])


    hbins=round((len(query_hits)/4))
    figure_to_png = plt.figure(figsize=[10,6])
    hist_to_png = plt.hist(query_hits,bins=hbins)
    hist_to_png = plt.axvline(sum(query_hits)/len(query_hits),color='k',linestyle='dashed',linewidth=2, label='average '+str(round(sum(query_hits)/len(query_hits))))
    hist_to_png = plt.axvline(min(query_hits),color='r',linestyle='dashed',linewidth=2, label='minimum '+str(min(query_hits)))
    hist_to_png = plt.axvline(max(query_hits),color='green',linestyle='dashed',linewidth=2, label='maximum '+str(max(query_hits)))
    #hist_to_png = plt.xticks(np.arange(0, max(query_hits), step=50))
    hist_to_png = plt.xlabel("amount of hits")
    hist_to_png = plt.ylabel("amount of query sequences")
    plt.savefig(snakemake.params['amount_hits_png'])
    plt.savefig(snakemake.output[2])

    figure_to_html = plt.figure(figsize=[10, 6])
    label = ['average ' + str(
        round(sum(query_hits) / len(query_hits))), 'minimum ' + str(
        min(query_hits)), 'maximum ' + str(max(query_hits))]
    hist_to_html = plt.hist(query_hits, bins=hbins, fc='green', label=label)
    hist_to_html = plt.grid(color='grey', linestyle='dashed')
    hist_to_html = plt.legend()
    hist_to_html = plt.xlabel("amount of hits")
    hist_to_html = plt.ylabel("amount of query sequences")
    mpld3.save_html(figure_to_html, snakemake.output[3])