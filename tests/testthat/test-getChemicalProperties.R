for (rampDB in test_databases) {
  test_that("chem props returns correctly formatted output for metabolites of interest, ChemicalProperties", {

    mets = c('hmdb:HMDB0000056',
             'hmdb:HMDB0000439'
    )
    chemProps <- getChemicalProperties(db = rampDB, mets = mets, propertyList = c('iso_smiles'))

    cn <- c('chem_source_id', 'iso_smiles')
    data <-
      c(
        'hmdb:HMDB0000056',
        'NCCC(O)=O',
        'hmdb:HMDB0000439',
        'OC(=O)CNC(=O)C1=CC=CO1'
      )
    df <- data.frame(matrix(
      data = data,
      ncol = 2,
      byrow = TRUE
    ))

    colnames(df) <- cn

    result <- chemProps$chem_props
    result <- result[, -2]
    #df <- df[, -2]
    df <- as.data.frame(df)
    result <- as.data.frame(result)
    print(result)
    print(df)

    #df <- df[order(df$chem_source_id), ]
    df <- df[order(as.vector(df$chem_source_id)),]
    result <- df[order(df$chem_source_id),]


    #print(df)
    print(table(result == df))
    res <- all.equal(result, df)
    #print(res)
    expect_equal(res, TRUE)
  })

}
