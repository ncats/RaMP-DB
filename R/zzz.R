
.onLoad <- function(libname, pkgname) {
}

.onAttach <- function(libname, pkgname) {

  print(paste0("Loading Ramp version: ", packageDescription("RaMP")$Version))
  print("Checking for existing BiocFileCache entry for the RaMP Database.")

  bfc <- BiocFileCache::BiocFileCache(cache = BiocFileCache::getBFCOption("CACHE"), ask=F)
  cacheInfo <- BiocFileCache::bfcinfo()
  cacheInfo <- cacheInfo[grepl("RaMP", cacheInfo$rname),]

  if(nrow(cacheInfo) < 1) {
    print("")
    print("RaMP Database is not in file cache. Performing one-time download.")
    url = packageDescription("RaMP")$Config_ramp_db_url
    print("One time retrieval of RaMP Database Cache. This will take about 1 minute to download and unzip.")
    path <- BiocFileCache::bfcadd(bfc, url, fname='exact')
    cid <- names(path)
    R.utils::gunzip(path, remove=F)
    newpath <- gsub(".gz", "", path)
    bfcremove(bfc, cid)
    bfcadd(bfc, newpath, fname='exact')
  } else {
    print(paste0("RaMP DB found in BiocFileCache, SQLite File: ",cacheInfo$rpath[1]))
  }
}



