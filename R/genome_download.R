#' Download genome metadata and sequence from NCBI
#'
#' @param genome_id The accession number of the genome (e.g., "ON400500.1")
#'
#' @return A data.frame containing genome metadata and sequence
#'
#' @importFrom magrittr %>%
#' @importFrom purrr set_names
#' @export
#'
#' @examples
#' \dontrun{
#' Download a single genome
#' fgenome_download("ON400500.1")
#' }
fgenome_download <- function(genome_id = "ON400500.1") {
  Sys.sleep(0.1)
  # config args
  args <- list(db = "nuccore", id = genome_id, retmode = "xml",
               api_key =  Sys.getenv("key_ncbi"))

  # retrieve data from ncbi
  con <- crul::HttpClient$new("https://eutils.ncbi.nlm.nih.gov")
  w <- con$get(path = "entrez/eutils/efetch.fcgi",
               query = Filter(f = Negate(is.null), x = args))


  # parse function
  xml_helper <- function(y, string) {
    xml2::xml_text(xml2::xml_find_first(y, string))
  }

  xml_helpers <- function(y, string) {
    xml2::xml_text(xml2::xml_find_all(y, string))
  }

  # parse content
  xml <- xml2::read_xml(w$parse("UTF-8"))

  w$raise_for_status()
  df <- lapply(xml2::xml_children(xml), function(z) {
    gitmp <- xml2::xml_text(xml2::xml_find_all(z, './GBSeq_other-seqids//GBSeqid'))
    gi <- strsplit(gitmp[length(gitmp)], "\\|")[[1]][2]
    acc <- xml_helper(z, './GBSeq_accession-version')
    def <- xml_helper(z, './GBSeq_definition')
    seq <- xml_helper(z, './GBSeq_sequence')
    seqlen <- xml_helper(z, './GBSeq_length')
    tax <- xml_helper(z, "./GBSeq_organism")
    taxonomy <- xml_helper(z, "./GBSeq_taxonomy")
    date <- xml_helper(z, "./GBSeq_create-date")
    date.sample <- xml_helper(z, './/GBQualifier[GBQualifier_name = "collection_date"]/GBQualifier_value')
    host <- xml_helper(z, './/GBQualifier[GBQualifier_name = "host"]/GBQualifier_value')
    iso.source <- xml_helper(z, './/GBQualifier[GBQualifier_name = "isolation_source"]/GBQualifier_value')
    country <- xml_helper(z, './/GBQualifier[GBQualifier_name = "country"]/GBQualifier_value')
    first.author <- xml_helper(z, './/GBReference[GBReference_reference = "1"]/GBReference_authors/GBAuthor')
    paper.title <- xml_helper(z, './/GBReference[GBReference_reference = "1"]/GBReference_title')
    journal <- xml_helper(z, './/GBReference[GBReference_reference = "1"]/GBReference_journal')

    # create data. frame
    data.frame(gi_no = gi,
               acc_no = acc,
               gene_desc = def,
               host = host,
               isolation_source = iso.source,
               country = country,
               date_sampled = date.sample,
               date_uploaded = date,
               paper_title = paper.title,
               journal = journal,
               first_author = first.author,
               length = seqlen,
               sequence = seq,
               stringsAsFactors = FALSE)
  }) %>%
    dplyr::bind_rows()

  return(df)

}

#' Download feature table for a genome from NCBI
#'
#' @param genome_id The accession number of the genome (e.g., "ON400500.1")
#'
#' @return A data.frame containing feature table information
#' @export
#'
#' @examples
#' \dontrun{
#' ft_download("ON400500.1")
#' }
ft_download <- function(genome_id = "ON400500.1") {
  Sys.sleep(0.5)
  # config args
  print(glue::glue("{genome_id}......"))

  args_ft <- list(db = "nuccore", id = genome_id, rettype = "ft",
                  api_key =  Sys.getenv("key_ncbi"))
  # retrieve data from ncbi
  con <- crul::HttpClient$new("https://eutils.ncbi.nlm.nih.gov")

  w_ft <- con$get(path = "entrez/eutils/efetch.fcgi",
                  query = Filter(f = Negate(is.null), x = args_ft))

  xml_ft <- readr::read_tsv(w_ft$parse("UTF-8"), skip = 1, col_names = FALSE, show_col_types = FALSE)
  if (nrow(xml_ft) > 0) {
    ft <- xml_ft %>%
      set_names(c("start", "end", "desc")) %>%
      tidyr::fill(start, .direction = "down") %>%
      tidyr::fill(end, .direction = "down") %>%
      tidyr::separate(col = desc, into = c("type", "name"), sep = "\\t", extra = "merge", fill = "right") %>%
      tidyr::drop_na()

    if (sum(stringr::str_detect(ft$type, "gene|protein_id")) > 0) {
      if (sum(stringr::str_detect(ft$type, "gene")) == 0){
        ft_w <- ft %>%
          dplyr::filter(stringr::str_detect(type,"gene|protein_id")) %>%
          dplyr::mutate(name = stringr::str_remove_all(name, pattern = "[gb|]")) %>%
          dplyr::mutate(start = as.numeric(start),
                        end = as.numeric(end)) %>%
          dplyr::mutate(cdsLen = abs(end - start) + 1) %>%
          dplyr::group_by(start, end) %>%
          dplyr::mutate(minpos = min(start, end)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(acc_no = genome_id) %>%
          dplyr::arrange(minpos, type)  %>%
          tidyr::pivot_wider(names_from = type, values_from = name) %>%
          dplyr::mutate(gene = "NA") %>%
          dplyr::select(acc_no, gene, protein_id, start, end, cdsLen, minpos) %>%
          tidyr::unnest(gene) %>%
          tidyr::unnest(protein_id)
        print(glue::glue("finished! But gene_id is unknown\n\n"))
        return(ft_w)
      } else if (sum(stringr::str_detect(ft$type, "protein_id")) == 0) {
        ft_w <- ft %>%
          dplyr::filter(stringr::str_detect(type,"gene|protein_id")) %>%
          dplyr::mutate(name = stringr::str_remove_all(name, pattern = "[gb|]")) %>%
          dplyr::mutate(start = as.numeric(start),
                        end = as.numeric(end)) %>%
          dplyr::mutate(cdsLen = abs(end - start) + 1) %>%
          dplyr::group_by(start, end) %>%
          dplyr::mutate(minpos = min(start, end)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(acc_no = genome_id) %>%
          dplyr::arrange(minpos, type)  %>%
          tidyr::pivot_wider(names_from = type, values_from = name) %>%
          dplyr::mutate(protein_id = "NA") %>%
          dplyr::select(acc_no, gene, protein_id, start, end, cdsLen, minpos) %>%
          tidyr::unnest(gene) %>%
          tidyr::unnest(protein_id)
        print(glue::glue("finished! But protein_id is unknown\n\n"))
        return(ft_w)
      } else {
        ft_w <- ft %>%
          dplyr::filter(stringr::str_detect(type,"gene|protein_id")) %>%
          dplyr::mutate(name = stringr::str_remove_all(name, pattern = "[gb|]")) %>%
          dplyr::mutate(start = as.numeric(start),
                        end = as.numeric(end)) %>%
          dplyr::mutate(cdsLen = abs(end - start) + 1) %>%
          dplyr::group_by(start, end) %>%
          dplyr::mutate(minpos = min(start, end)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(acc_no = genome_id) %>%
          dplyr::arrange(minpos, type)  %>%
          tidyr::pivot_wider(names_from = type, values_from = name) %>%
          dplyr::select(acc_no, gene, protein_id, start, end, cdsLen, minpos) %>%
          tidyr::unnest(gene) %>%
          tidyr::unnest(protein_id)
        message(glue::glue("Successfully downloaded feature table for {genome_id}"))
        return(ft_w)
      }
    } else {
      message(glue::glue("No gene_id or protein_id found for {genome_id}!"))
      return(NULL)
      }
  }else {
    message(glue::glue("Failed to download feature table for {genome_id}!"))
    return(NULL)  }
}
