library(tidyverse)
read_salmon <- function(path, which_return='tibble',quant_type='counts', qfiles = '', normalize_counts=T){
    if(qfiles[1] == ''){
        qfiles <- list.files(path,pattern='quant.sf', recursive=T, full.names = T)
    }
    name_idx <- str_split(qfiles[1], '/')[[1]] %>% grep('quant.sf', .) %>% {. -1}
    names <- str_split(qfiles,'/') %>% sapply(function(x) x[name_idx])
    if(normalize_counts){
        txi <- tximport::tximport(files = qfiles, type = 'salmon', txOut = T, countsFromAbundance = 'lengthScaledTPM')
    }else{
        
        txi <- tximport::tximport(files = qfiles, type = 'salmon', txOut = T)
    }
    colnames(txi$counts) <- names
    colnames(txi$abundance) <- names
    if(which_return == 'tibble'){
        
	res <- txi[[quant_type]] %>% as.data.frame %>% dplyr::mutate(transcript_id =rownames(.)) %>% dplyr::select(transcript_id, everything())
        return(res)
    } else if (which_return == 'txi'){
            return(txi)
    } else{
        message('bad input')
    }
}
