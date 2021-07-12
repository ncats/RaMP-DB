test_that("multiplication works", {
  loc <- getwd()
  print(loc)
  dbpass <-
    read.csv(
      file = '../../dbprops.txt',
      sep = ",",
      header = FALSE,
      quote = '"',
      row.names = 1
    )




  mets = c('hmdb:HMDB0000056',
           'hmdb:HMDB0000439',
           'hmdb:HMDB0000479')

  chemProps <-
    getChemicalProperties(
      mets,
      propertyList = c('smiles'),
      conpass = dbpass['conpass', 1],
      dbname = dbpass['dbname', 1],
      host = dbpass['host', 1],
      username = dbpass['username', 1]
    )
  print(head(chemProps$chem_props))


  cn <- c('chem_source_id', 'ramp_id', 'iso_smiles')
  data <-
    c(
      'hmdb:HMDB0000056',
      'RAMP_C_000088275',
      'NCCC(O)=O',
      'hmdb:HMDB0000479',
      'RAMP_C_000055801',
      'CN1C=NC=C1C[C@H](N)C(O)=O',
      'hmdb:HMDB0000439',
      'RAMP_C_000043583',
      'OC(=O)CNC(=O)C1=CC=CO1'
    )
  df <- data.frame(matrix(
    data = data,
    ncol = 3,
    byrow = TRUE
  ))
  colnames(df) <- cn
  print(df)
  print(table(chemProps$chem_props == df))
  res <- all.equal(chemProps$chem_props, df)
  print(res)



  expect_equal(res, TRUE)
})
