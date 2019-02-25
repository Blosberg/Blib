# ==================================================================
# --- get arbitrary sequences from an established reference

get_refgen_seqs <-  function(  refgen     = stop("seq  must be provided"),
                               ROI_GR     = stop("Region of interest must be provided in GRanges format"),
                               lead       = 0,
                               trail      = 0,
                               RNAstrand  = FALSE
)
{
  if ( as.character(strand(ROI_GR) ) == "+")
  {

    lo_end   <- start(ROI_GR) - lead;
    hi_end   <- end(ROI_GR)   + trail;
    flipseq = FALSE;

  }else if( as.character(strand(ROI_GR) )  == "-") {

    lo_end   <- start(ROI_GR) - trail;
    hi_end   <- end(ROI_GR)   + lead;    # it will be flipped later.
    flipseq = TRUE;

  } else {
    stop("undefined strand in get_refgen_seqs");
  }

  result = eval (
                parse(
                     text=paste0(
                                "refgen$",seqnames(ROI_GR),"[",as.character(lo_end),":",as.character(hi_end),"]"
                                )
                    )
                )

  # ------------------------------
  if ( flipseq )
  { # --- if the read is from the - strand, then reverse and compliment it
    result = chartr("ATGC","TACG", stringi::stri_reverse( result ) )
  }

  if( RNAstrand )
  {
    result = sub( "T", "U", result )
  }

  return( result )
}

# ==================================================================
# --- Extend sequence provided in pore data, from input reference.

expand_range <- function ( GRin       = stop("GRmod must be provided"),
                           refGenome  = stop("GRcontrol must be provided"),
                                    k = 5 )
{

  GR_justlocs = GRanges( seqnames = seqnames(GRin),
                         ranges   = IRanges(start = start(GRin), ))



  if ( k <= 0)
  {stop("k must be a positive integer")}
  # --------

  if ( as.character(strand(GRin) ) == "+")
  {
    lo_end   <- start(GRin) ;
    hi_end   <- start(GRin) -1 + k;
    flipseq = FALSE;

  }else if( as.character(strand(ROI_GR) )  == "-") {

    lo_end   <- start(ROI_GR) - trail;
    hi_end   <- end(ROI_GR)   + lead;    # it will be flipped later.
    flipseq = TRUE;

  }

  result = eval (
    parse(
      text=paste0(
        "refGenome$",seqnames(GRin),"[",as.character(lo_end),":",as.character(hi_end),"]"
      )
    )
  )

 return(result)
}


# ==================================================================
# --- Go through a set of read sequences and reverse the order
# --- (and kmer) of only the ones aligned to the reverse strand

strand_align <- function( GRin = stop("GRobj must be provided"),
                          k = 6)
{ # k is the length of the kmer sequence (minION devices spit out k=5, but we can extend that further.)
  temp = GRin[strand(GRin)=="-"]

  GRin[ strand(GRin)=="-"]  <- strand_reversal( temp )

  return(GRin)
}



# ==================================================================
# --- reverse the strand of a given read

strand_reversal  <- function( range_in = stop("range to be flipped must be provided"),
                              k = 5 )
{
  temp                <- range_in
  temp$reference_kmer <- reverse( chartr("ATGC","TACG", range_in$reference_kmer ) )

  if( !is.null( temp$modposition) ) # i.e. check if the metadata column exists.
    {
    temp$modposition    <- (k+1) -range_in$modposition
    }

  return(temp)
}
