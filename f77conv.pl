#!/usr/bin/perl

my $line1;
my $line2;

$line1=<>;
chomp $line1;

while($line2=<>) {
   chomp $line2;
   if ( ($line2 !~ /^(!|[cC])/) && (length ($line2) >6) && substr($line2,5,1) ne ' ') {
      substr($line2,5,1)=' ';
      # but we need to look for a comment first
       if ($line1 =~ /^[^!](.*?)!/   ) {
         $line1 =~ s/!/& !/;
       }else{
         chomp $line1;
         $line1 .= " &";
      }
          
   }
   print "$line1\n";
   $line1=$line2;
}
print "$line2\n";

