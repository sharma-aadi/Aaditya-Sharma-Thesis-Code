## 1. Construct SummarizedExperiement
pe_dda = make_pe_parse(peptide_table_dda, 
                       columns = 5:43,   # columns is the 'Intensity' colunmns
                       mode = "delim", sep = "_",
                       remove_prefix = FALSE,
                       remove_suffix = FALSE
)

pe_dia = make_pe_parse(peptide_table_dia, 
                       columns = 5:43,   # columns is the 'Intensity' colunmns
                       mode = "delim", sep = "_",
                       remove_prefix = FALSE,
                       remove_suffix = FALSE
)

## 2. Filter
pe_dda = filter_pe(pe_dda, thr = 2, fraction = 0.1)
pe_dia = filter_pe(pe_dia, thr = 2, fraction = 0.1)

## 3. Impute

# No imputation
pe_dda_IN = pe_dda # peptideRaw will be the assay
pe_dia_IN = pe_dia # peptideRaw will be the assay

# QRILC
pe_dda_IQ = impute_pe(pe_dda, fun = "QRILC", name = "peptideImp")
pe_dia_IQ = impute_pe(pe_dia, fun = "QRILC", name = "peptideImp")

# MinProb
pe_dda_IM = impute_pe(pe_dda, fun = "MinProb", q = 0.05, name = "peptideImp")
pe_dia_IM = impute_pe(pe_dia, fun = "MinProb", q = 0.05, name = "peptideImp")

## 4. Normalize
pe_dda_IN = normalize_pe(pe_dda_IN, method = "diff.median", i = "peptideRaw") # peptideRaw will be the assay
pe_dia_IN = normalize_pe(pe_dia_IN, method = "diff.median", i = "peptideRaw") # peptideRaw will be the assay

pe_dda_IQ = normalize_pe(pe_dda_IQ, method = "diff.median", i = "peptideImp")
pe_dia_IQ = normalize_pe(pe_dia_IQ, method = "diff.median", i = "peptideImp")

pe_dda_IM = normalize_pe(pe_dda_IM, method = "diff.median", i = "peptideImp")
pe_dia_IM = normalize_pe(pe_dia_IM, method = "diff.median", i = "peptideImp")

## 5. Aggregate peptides quantity to proteins

### Robust Summary
begin_time = Sys.time()

pe_dda_IN_RS <- aggregate_pe(pe_dda_IN, fcol = "Leading.razor.protein", reserve = "Gene.names",
                          aggregate_Peptide_Type = "Unique", aggrefun = "RobustSummary")
pe_dia_IN_RS <- aggregate_pe(pe_dia_IN, fcol = "Leading.razor.protein", reserve = "Gene.names",
                          aggregate_Peptide_Type = "Unique", aggrefun = "RobustSummary")

pe_dda_IQ_RS <- aggregate_pe(pe_dda_IQ, fcol = "Leading.razor.protein", reserve = "Gene.names",
                          aggregate_Peptide_Type = "Unique", aggrefun = "RobustSummary")

pe_dia_IQ_RS <- aggregate_pe(pe_dia_IQ, fcol = "Leading.razor.protein", reserve = "Gene.names",
                          aggregate_Peptide_Type = "Unique", aggrefun = "RobustSummary")

pe_dda_IM_RS <- aggregate_pe(pe_dda_IM, fcol = "Leading.razor.protein", reserve = "Gene.names",
                          aggregate_Peptide_Type = "Unique", aggrefun = "RobustSummary")
pe_dia_IM_RS <- aggregate_pe(pe_dia_IM, fcol = "Leading.razor.protein", reserve = "Gene.names",
                             aggregate_Peptide_Type = "Unique", aggrefun = "RobustSummary")

print(Sys.time() - begin_time) ## Required 16.5 minutes

### Median polish
begin_time = Sys.time()

pe_dda_IN_MP <- aggregate_pe(pe_dda_IN, fcol = "Leading.razor.protein", reserve = "Gene.names",
                             aggregate_Peptide_Type = "Unique", aggrefun = "medianPolish")
pe_dia_IN_MP <- aggregate_pe(pe_dia_IN, fcol = "Leading.razor.protein", reserve = "Gene.names",
                             aggregate_Peptide_Type = "Unique", aggrefun = "medianPolish")

pe_dda_IQ_MP <- aggregate_pe(pe_dda_IQ, fcol = "Leading.razor.protein", reserve = "Gene.names",
                             aggregate_Peptide_Type = "Unique", aggrefun = "medianPolish")

pe_dia_IQ_MP <- aggregate_pe(pe_dia_IQ, fcol = "Leading.razor.protein", reserve = "Gene.names",
                             aggregate_Peptide_Type = "Unique", aggrefun = "medianPolish")

pe_dda_IM_MP <- aggregate_pe(pe_dda_IM, fcol = "Leading.razor.protein", reserve = "Gene.names",
                             aggregate_Peptide_Type = "Unique", aggrefun = "medianPolish")
pe_dia_IM_MP <- aggregate_pe(pe_dia_IM, fcol = "Leading.razor.protein", reserve = "Gene.names",
                             aggregate_Peptide_Type = "Unique", aggrefun = "medianPolish")

print(Sys.time() - begin_time) ## Required 1.6 minutes

## 6. Transform a SummarizedExperiement of protein quantities.
pe_dda_IN_RS_se <- pe2se(pe_dda_IN_RS, names = "Leading.razor.protein", ids = "smallestProteingroups")
pe_dia_IN_RS_se <- pe2se(pe_dia_IN_RS, names = "Leading.razor.protein", ids = "smallestProteingroups")

pe_dda_IQ_RS_se <- pe2se(pe_dda_IQ_RS, names = "Leading.razor.protein", ids = "smallestProteingroups")
pe_dia_IQ_RS_se <- pe2se(pe_dia_IQ_RS, names = "Leading.razor.protein", ids = "smallestProteingroups")

pe_dda_IM_RS_se <- pe2se(pe_dda_IM_RS, names = "Leading.razor.protein", ids = "smallestProteingroups")
pe_dia_IM_RS_se <- pe2se(pe_dia_IM_RS, names = "Leading.razor.protein", ids = "smallestProteingroups")

pe_dda_IN_MP_se <- pe2se(pe_dda_IN_MP, names = "Leading.razor.protein", ids = "smallestProteingroups")
pe_dia_IN_MP_se <- pe2se(pe_dia_IN_MP, names = "Leading.razor.protein", ids = "smallestProteingroups")

pe_dda_IQ_MP_se <- pe2se(pe_dda_IQ_MP, names = "Leading.razor.protein", ids = "smallestProteingroups")
pe_dia_IQ_MP_se <- pe2se(pe_dia_IQ_MP, names = "Leading.razor.protein", ids = "smallestProteingroups")

pe_dda_IM_MP_se <- pe2se(pe_dda_IM_RS, names = "Leading.razor.protein", ids = "smallestProteingroups")
pe_dia_IM_MP_se <- pe2se(pe_dia_IM_MP, names = "Leading.razor.protein", ids = "smallestProteingroups")

#######################################################################################################################
## Save
#######################################################################################################################

# save.image("LFQProBe_DEA_datasets.RData")
