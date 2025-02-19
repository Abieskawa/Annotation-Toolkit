use strict;
use warnings;
use LWP::UserAgent;
use LWP::Simple;
use HTTP::Date;

# The script was copied from https://www.uniprot.org/help/api_downloading
# Taxonomy identifier of top node for query, e.g. 2 for Bacteria, 2157 for Archea, etc.
# (see https://www.uniprot.org/taxonomy)
my $top_node = $ARGV[0];

my $agent = LWP::UserAgent->new;

# Get TSV of all reference proteomes of organisms below the given taxonomy node.
my $query_list = "https://rest.uniprot.org/proteomes/stream?&query=reference:true+taxonomy_id:$top_node&fields=upid,lineage,organism_id&format=tsv";

my $response_list = $agent->get($query_list);
if ( $response_list->is_success ) {
    my $release = $response_list->header('x-uniprot-release');
    my $date    = $response_list->header('x-uniprot-release-date');
    print "Fetching FASTAs from UniProt release $release ($date)\n";
}
elsif ( $response_list->code == HTTP::Status::RC_NOT_MODIFIED ) {
    print "Data for taxon $top_node is up-to-date.\n";
}
else {
    die 'Failed, got '
      . $response_list->status_line . ' for '
      . $response_list->request->uri . "\n";
}

# For each proteome, fetch compressed FASTA from FTP.
my @lines = split( /\n/, $response_list->content );
# Skip the TSV header by starting from 1.
for my $index ( 1 .. $#lines ) {
    my @line                     = split( /\t/, $lines[$index] );
    my $upid                     = $line[0];
    my @taxonomic_lineage_column = split( /,\s/, $line[1] );
    my $domain                   = $taxonomic_lineage_column[1]; # first column is "cellular organisms", second is kingdom
    my $organism_id              = $line[2];
    my $url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/$domain/$upid/$upid\_$organism_id.fasta.gz";
    my $file = "$upid.fasta.gz";
    my $response_ftp = getstore( $url, $file );

    if ( is_success($response_ftp) ) {
        print "Fetched $url\n";
    }
    if ( is_error($response_ftp) ) {
        die "getstore of <$url> failed with $response_ftp";
    }
}