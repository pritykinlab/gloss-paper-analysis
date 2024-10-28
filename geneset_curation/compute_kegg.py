def get_kegg_pathways_global():
    # Initialize KEGG service
    kegg = KEGG()
    # Specify the organism for mouse
    organism = "mmu"
    # Get the list of all mouse pathways
    pathways = kegg.list("pathway", organism)
    
    # Extract pathway IDs
    pathway_ids = [line.split()[0] for line in pathways.splitlines()]

    # Initialize an empty dictionary to store the pathways
    pathway_dict = {}

    # Loop through each pathway ID and store the data in the dictionary
    for i, pathway_id in enumerate(pathway_ids):
        if i % 50 == 0:
            print(i)
        data = kegg.get(pathway_id)
        pathway_dict[pathway_id] = data
    
    pathway_genes_dict = {}

    # Loop through each pathway ID and extract the genes
    for pathway_id in pathway_ids:
        data = pathway_dict[pathway_id]
    
        # Parse the genes from the pathway data (look for entries starting with 'GENE')
        genes = []
        recording_genes = False
        for line in data.splitlines():
            if line.startswith("GENE"):
                recording_genes = True
                gene_info = line.split()[2][:-1]  # Split and take the gene identifier
                genes.append(gene_info)
            elif recording_genes:
                if line.startswith(" "):
                    gene_info = line.split()[1][:-1]  # Split and take the gene identifier
                    genes.append(gene_info)
                else:
                    break

        # Store the set of genes for this pathway in the dictionary
        pathway_genes_dict[pathway_id] = genes
        
    pathway_id_names = [line.split('\t')[1].split(' - Mus musculus')[0] for line in pathways.splitlines()]

    newdict = {}
    for i, pathway_id in enumerate(pathway_genes_dict):
        newdict[pathway_id_names[i]] = pathway_genes_dict[pathway_id]

    kegg_gene_sets = newdict
        
    empty = []
    for key in kegg_gene_sets:
        if not kegg_gene_sets[key]:
            empty.append(key)
        
    path_dict = kegg_gene_sets
    return path_dict

