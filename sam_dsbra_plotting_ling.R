#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))

# read in arguments
option_list <- list(
    make_option('--input', help='input file'),
    make_option('--mut_type', action='store_true', default=FALSE,
                help="plot mutation_event_frequency_by_type.csv,\
                      plot Mutation Event Frequency by Type"),
    make_option("--del_len", action='store_true', default=FALSE,
                help="plot deletion_lens.csv,\
                plot Frequency of Deletions by Length"),
    make_option("--del_seq", action='store_true', default=FALSE,
                help="plot sequences_with_deletion_events.csv\
                plot Sequences with Deletion Event"),
    make_option("--min_count", default=0,
                help="determine how the minimum counts to include in\
                      Sequences with Deletion/Insertion/Mutation Event"),
    make_option("--ins_len", action='store_true', default=FALSE,
                help="plot insertion_lens.csv\
                plot Frequency of Insertions by Length"),
    make_option("--ins_seq", action='store_true', default=FALSE,
                help="plot sequences_with_insertion_events.csv\
                plot Sequences with Insertion Event"),
    make_option("--repair_seq", action='store_true', default=FALSE,
                help="plot sequences_with_mutation_events.csv\
                plot Sequences with Mutation Event"),
    make_option("--aligned_mutations", action='store_true', default=FALSE,
                help="plot aligned_mutation_events.csv\
                plot mutation events aligned with the reference sequence"),
    make_option("--ref_bottom", default=0, type="integer",
                help="the start reference base to show aligned mutation events"),
    make_option("--ref_top", type="integer",
                help="the end reference base to show aligned mutation events"),
    make_option("--break_index", type="integer",
                help="the index of break site")
    )
opt <- parse_args(OptionParser(option_list=option_list))

# read input dataframe
df <- read.csv(opt$input)
basename <- tools::file_path_sans_ext(opt$input)

# output name
outname <- sprintf('%s.png', basename)
print(sprintf('Writing to %s', outname))

min_count <- opt$min_count

# Mutation Event Frequency by Type, pie
#if (opt$mut_type){
#    png(outname)
#    dft <- t(data.frame(count=df$count, row.names = df$mutation_type))
#    lbls <- paste(colnames(dft), "\n", df$count)
#    pie(dft, labels = lbls,
#        main="Mutation Event Frequency by Type")
#    dev.off()
#}

# Mutation Event Frequency by Type, bar
if (opt$mut_type){
#    df$mutation_type <- factor(df$mutation_type, levels = df[order(df$count), ]$mutation_type)
    # order by deletion insertion compound WT
     df$mutation_type <- factor(df$mutation_type, levels = c("Deletion", "Insertion", "Compound", "WT"))
    p <- ggplot(data=df, aes(x=mutation_type, y=count)) +
         geom_bar(stat="identity", width=.5) +
         labs(x='Mutation Type', y='Count', title='Mutation Event Frequency by Type') +
         geom_text(aes(label=count), position=position_dodge(width=0.9), vjust=-0.3, size=3) +
         theme_bw() +
         theme(aspect.ratio = 0.6)
    ggsave(filename=outname, plot=p)
}

# Frequency of Deletions by Length
if (opt$del_len){
    p <- ggplot(df, aes(x=deletion_length)) +
        geom_histogram(colour="black") +
        labs(x='Deletion Length (bp)', y='Count', title='Frequency of Deletions by Length') +
        theme_bw() +
        theme(aspect.ratio = 0.6)
    ggsave(filename=outname, plot=p)
}

# Sequences with deletion events
if (opt$del_seq){
    df$sequence <- factor(df$sequence, levels = df[order(df$count),]$sequence)
    df <- subset(df, count > min_count)
    p <- ggplot(data=df[order(-df$count),], aes(x=sequence, y=count)) +
    geom_bar(stat="identity", width=.5) +
    labs(x='Sequences', y='Count', title='Sequences With Deletion Events') +
    geom_text(aes(label=count), position=position_dodge(width=0.9), hjust=-0.3, size=4*20/nrow(df)) +
    theme_bw() +
    theme(aspect.ratio = 1 * nrow(df) / 20,
          axis.text=element_text(size=14*20/nrow(df))) +
    coord_flip(ylim = c(0, max(df$count) * 1.2))
    ggsave(filename=outname, plot=p)
}


# Frequency of Insertion by Length
if (opt$ins_len){
    p <- ggplot(df, aes(x=insertion_length)) +
        geom_histogram(colour="black") +
        labs(x='Insertion Length (bp)', y='Count', title='Frequency of Insertions by Length') +
        theme_bw() +
        theme(aspect.ratio = 0.6)
    ggsave(filename=outname, plot=p)
}

# Sequences with insertion events
if (opt$ins_seq){
    df$sequence <- factor(df$sequence, levels = df[order(df$count),]$sequence)
    df <- subset(df, count > min_count)
    p <- ggplot(data=df[order(-df$count),], aes(x=sequence, y=count)) +
    geom_bar(stat="identity", width=.5) +
    labs(x='Sequences', y='Count', title='Sequences With Insertion Events') +
    geom_text(aes(label=count), position=position_dodge(width=0.9), hjust=-0.3, size=3) +
    theme_bw() +
    theme(aspect.ratio = 0.5 * nrow(df) / 20) +
    coord_flip(ylim = c(0, max(df$count) * 1.2))
    tryCatch(ggsave(filename=outname, plot=p), error= function(err){print('insetion seq fails')})
}

# Sequences with mutation events
if (opt$repair_seq){
    df$sequence <- factor(df$sequence, levels = df[order(df$count),]$sequence)
    df <- subset(df, count > min_count)
    p <- ggplot(data=df[order(-df$count),], aes(x=sequence, y=count)) +
    geom_bar(stat="identity", width=.5) +
    labs(x='Sequences', y='Count', title='Repair Patterns') +
    geom_text(aes(label=count), position=position_dodge(width=0.9), hjust=-0.3, size=4*20/nrow(df)) +
    theme_bw() +
    theme(aspect.ratio = 1 * nrow(df) / 20,
          axis.text=element_text(size=14*20/nrow(df))) +
    coord_flip(ylim = c(0, max(df$count) * 1.2))
    tryCatch(ggsave(filename=outname, plot=p), error= function(err){print('repair patterns fails')})
}

# Aligned mutation events
if (opt$aligned_mutations){
    ref_seq_range <- c(opt$ref_bottom, opt$ref_top)
    break_index <- opt$break_index
    mut_size_font <- 1.5
    mut_size_hjust <- -3
    aspect_ratio <- 1
    width_scale <- nrow(df) * 1/2000

    # prepare table
    df <- df[, 2:8]
    # because if mut_start = 85, then it will cover [85,86].
    # If break site is 84, then it is marked at [83, 84] | [84,85].
    # so we want to mark the mutation happened at 85, as [84, 85]
    df$mut_start <- df$mut_start - break_index - 1
    df$mut_end <- df$mut_end - break_index - 1

    x_start <- ref_seq_range[1] - break_index
    x_end <- ref_seq_range[2] - break_index
    x_interval <- 5
    # make sure 0bp is shown on the plot
    while ((0 - x_start) %% 5){
        x_start <- x_start - 1
    }

    df2 <- data.frame(idx_start = df$idx_start,
                      idx_end = df$idx_end,
                      mut_start = rep(ref_seq_range[1], nrow(df)) - break_index,
                      mut_end = rep(ref_seq_range[2], nrow(df)) - break_index,
                      type = rep('matched', nrow(df)),
                      region_len = NA,
                      count = df$count)

    dft <- rbind(df2, df)
    p1 <- ggplot(dft, aes(xmin = mut_start, xmax = mut_end,
                          ymin = idx_start*width_scale,
                          ymax = idx_end*width_scale)) +
                 geom_rect(aes(fill = type)) +
                 geom_vline(xintercept = 0, color="#FF2E4C", linetype=3) +
                 #geom_text(aes(x = mut_start+mut_size_hjust,
                 #              y = (idx_end - idx_start)/2 + idx_start,
                 #          label=region_len),
                 #          size=mut_size_font) +
                 labs(x="bp from cut site", y="") +
                 guides(fill=guide_legend(title="Mutation Type")) +
                 scale_x_continuous(breaks=seq(x_start, x_end, x_interval)) +
                 scale_y_continuous(breaks=seq(0, 1.02*max(dft$idx_end), 1)) +  # cut size text margin 1.02 * max(y)
#                 annotate("text", label = "cut site", x = 0 ,
#                          y =  1.03*max(dft$idx_end), color = "#FF2E4C") +
                 scale_fill_manual(values = c("matched"="grey",
                                              "deletion"="black",
                                              "insertion"="#2E99B0",
                                              "mismatch"="#FCD77F")) +
                 theme(legend.position="bottom",
                       aspect.ratio =  aspect_ratio,
                       panel.background = element_blank(),
                       axis.text.x= element_text(size=8),
                       axis.line.x =element_line(color="black"),
                       axis.text.y= element_blank(),
                       axis.ticks.y=element_blank())
    ggsave(filename=outname, plot=p1)
}
