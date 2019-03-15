#' @importFrom dplyr bind_cols

eliminate <- function(MFs,by,direction){
  MFs %>%
    bind_cols(MFs %>% select(by = by)) %>%
    split(.$Feature) %>%
    map(~{
      d <- .
      direct <- get(direction)
      d %>%
        filter(by == direct(by))
    }) %>%
    bind_rows() %>%
    select(-by)
}