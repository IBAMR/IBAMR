# $Id: Simple.pm 2199 2007-03-17 01:08:55Z comdog $
package ConfigReader::Simple;
use strict;

use subs qw(_init_errors);
use vars qw($VERSION $AUTOLOAD %ERROR $ERROR $Warn $Die);

use Carp qw(croak carp);
use UNIVERSAL qw(isa);

$Die   = '';
$ERROR = '';
( $VERSION ) = 1.25;
#= sprintf "%d.%02d", q$Revision: 2199 $ =~ m/ (\d+) \. (\d+) /gx;
$Warn = 0;

our $DEBUG = 0;
my $Error = '';

sub SUCCESS() { 1 };
sub FAILURE() { 0 };

=head1 NAME

ConfigReader::Simple - Simple configuration file parser

=head1 SYNOPSIS

	use ConfigReader::Simple;

	# parse one file
	$config = ConfigReader::Simple->new("configrc", [qw(Foo Bar Baz Quux)]);

	# parse multiple files, in order
	$config = ConfigReader::Simple->new_multiple(
		Files => [ "global", "configrc" ], 
		Keys  => [qw(Foo Bar Baz Quux)]
		);

	my @directives = $config->directives;

	$config->get( "Foo" );

	if( $config->exists( "Bar" ) )
   		{
   		print "Bar was in the config file\n";
   		}
   		
   	# copy an object to play with it separately
   	my $clone = $config->clone;
   	
   	# only affects clone
   	$clone->set( "Foo", "Buster" );

	# save the config to a single file
	$clone->save( "configrc" )

	# save the config to a single file, but only with
	# certain directives
	$clone->save( "configrc" => [qw(Foo Bar)] )
	
	# save to multiple configuration files
	$clone->save( 
		"configrc" => [qw(Foo Bar)],
		"global"   => [qw(Baz Quux)],
		);

=head1 DESCRIPTION

C<ConfigReader::Simple> reads and parses simple configuration files. It is
designed to be smaller and simpler than the C<ConfigReader> module
and is more suited to simple configuration files. 

=head2 The configuration file format

The configuration file uses a line-oriented format, meaning
that the directives do not have containers.  The values can
be split across lines with a continuation character, but for
the most part everything ends up on the same line.

The first group of non-whitespace characters is the
"directive", or the name of the configuration item.  The
linear whitespace after that separates the directive from
the "value", which is the rest of the line, including any
other whitespace.

In this example, the directive is "Camel" and the value is
"Dromedary".

	Camel Dromedary
	
Optionally, you can use a equal sign to separate the directive
from the value.

	Camel=Dromedary
	
The equal sign can also have whitespace on either or both
sides.

	Camel = Dromedary
	Camel= Dromedary
	
In the next example, the directive is "Llama" and the value
is "Live from Peru"

	Llama Live from Peru
	
This is the same, to ConfigReader::Simple, as the following
which has more whitespace between the directive and the value.

	Llama     Live from Peru
	
You can also enclose the value in single or double quotes.

	Llama "Live from Peru"
	Llama 'Live from Peru'
	Llama='Live from Peru'

In some cases you may want to split the logical line across
two lines, perhaps to see it better in a terminal window.
For that, use a \ followed only by whitespace.  To split the
last entry across two lines, we use the \ at the end of the
line. These three entries are the same:

	Llama Live from Peru

	Llama Live from \
	Peru

	Llama Live \
	from \
	Peru

If a line is only whitespace, or the first whitespace character is 
a #, the Perl comment character, ConfigReader::Simple ignores the
line unless it is the continuation of the previous line.

=head2 Methods

=over 4

=item new ( FILENAME, DIRECTIVES )

Creates a ConfigReader::Simple object.

C<FILENAME> tells the instance where to look for the
configuration file. If FILENAME cannot be found, an error
message for the file is added to the %ERROR hash with the
FILENAME as a key, and a combined error message appears in
$ERROR.

C<DIRECTIVES> is an optional argument and is a reference to
an array. Each member of the array should contain one valid
directive. A directive is the name of a key that must occur
in the configuration file. If it is not found, the method
croaks. The directive list may contain all the keys in the
configuration file, a sub set of keys or no keys at all.

The C<new> method is really a wrapper around C<new_multiple>.

=cut

sub new 
	{
	my $class    = shift;
	my $filename = shift;
	my $keyref   = shift;
	
	$keyref = [] unless defined $keyref;
	
	my $self = $class->new_multiple( 
		Files => [ defined $filename ? $filename : () ],
		Keys  => $keyref );
			
	return $self;
	}

=item new_multiple( Files => ARRAY_REF, Keys => ARRAY_REF )

Create a configuration object from several files listed
in the anonymous array value for the C<Files> key.  The
module reads the files in the same order that they appear
in the array.  Later values override earlier ones.  This
allows you to specify global configurations which you 
may override with more specific ones:

	ConfigReader::Simple->new_multiple(
		Files => [ qw( /etc/config /usr/local/etc/config /home/usr/config ) ],
		);

This function croaks if the values are not array references.

If this method cannot read a file, an error message for that
file is added to the %ERROR hash with the filename as a key,
and a combined error message appears in $ERROR.  Processing
the list of filenames continues if a file cannot be found,
which may produced undesired results. You can disable this
feature by setting the $ConfigReader::Simple::Die variable
to a true value.

=cut

sub new_multiple
	{
	_init_errors();
	
	my $class    = shift;
	my %args     = @_;

	my $self = {};
	
	$args{'Keys'} = [] unless defined $args{'Keys'};
	
	croak( __PACKAGE__ . ': Strings argument must be a array reference')
		unless UNIVERSAL::isa( $args{'Files'}, 'ARRAY' );
	croak( __PACKAGE__ . ': Keys argument must be an array reference')
		unless UNIVERSAL::isa( $args{'Keys'}, 'ARRAY' );
		
	$self->{"filenames"} = $args{'Files'};
	$self->{"validkeys"} = $args{'Keys'};
	
	bless $self, $class;
	
	foreach my $file ( @{ $self->{"filenames"} } )
		{
		my $result = $self->parse( $file );
		croak $Error if( not $result and $Die );
		
		$ERROR{$file} = $Error unless $result;
		}
		
	$ERROR = join "\n", map { $ERROR{$_} } keys %ERROR;
	
	return $self;
	}

=item new_string( Strings => ARRAY_REF, Keys => ARRAY_REF )

Create a configuration object from several strings listed
in the anonymous array value for the C<Strings> key.  The
module reads the strings in the same order that they appear
in the array.  Later values override earlier ones.  This
allows you to specify global configurations which you 
may override with more specific ones:

	ConfigReader::Simple->new_strings(
		Strings => [ \$global, \$local ],
		);

This function croaks if the values are not array references.

=cut

sub new_string
	{
	_init_errors;
	
	my $class = shift;
	my %args  = @_;
	
	my $self = {};
	
	$args{'Keys'} = [] unless defined $args{'Keys'};

	croak( __PACKAGE__ . ': Strings argument must be a array reference')
		unless UNIVERSAL::isa( $args{'Strings'}, 'ARRAY' );
	croak( __PACKAGE__ . ': Keys argument must be an array reference')
		unless UNIVERSAL::isa( $args{'Keys'}, 'ARRAY' );

	bless $self, $class;

	$self->{"strings"} = $args{'Strings'};
	$self->{"validkeys"} = $args{'Keys'};
	
	foreach my $string ( @{ $self->{"strings"} } )
		{
		$self->parse_string( $string );
		}
		
	return $self;
	}
	
=item add_config_file( FILENAME )

Parse another configuration file and add its directives to the
current configuration object. Any directives already defined 
will be replaced with the new values found in FILENAME.

=cut

sub add_config_file
	{
	_init_errors;
	
	my $self     = shift;
	my $filename = shift;
	
	return unless ( -e $filename and -r _ );
	
	push @{ $self->{"filenames"} }, $filename
		if $self->parse( $filename );
	
	return 1;
	}

=item files

Return the list of configuration files associated with this 
object. The order of the return values is the order of parsing,
so the first value is the first file parsed (and subsequent files may
mask it).

=cut

sub files { @{ $_[0]->{"filenames"} } }

=item new_from_prototype( 

Create a clone object. This is the same thing as calling
clone().

=cut

sub new_from_prototype
	{
	_init_errors;

	my $self     = shift;
	
	my $clone = $self->clone;
	
	return $clone;
	}
	
sub AUTOLOAD
	{
	my $self = shift;

	my $method = $AUTOLOAD;

	$method =~ s/.*:://;

	$self->get( $method );
	} 

sub DESTROY 
	{	
	return 1;
	}

=item parse( FILENAME )

This does the actual work.

This is automatically called from C<new()>, although you can reparse
the configuration file by calling C<parse()> again.

=cut

sub parse 
	{
	my $self = shift;
	my $file = shift;
	
	$Error = '';
	
	unless( open CONFIG, $file )
		{
		$Error = "Could not open configuration file [$file]: $!";
		
		carp "Could not open configuration file [$file]: $!" if
			$Warn;
			
		return;
		}
		
	$self->{"file_fields"}{$file} = [];
	
	while( <CONFIG> )
		{
		if ( s/\\ \s* $//x )
			{
			$_ .= <CONFIG>;
			redo unless eof CONFIG;
			}
			
		chomp;
		next if /^\s*(#|$)/; 
		
		my ($key, $value) = &parse_line($_);
		carp "Key:  '$key'   Value:  '$value'\n" if $DEBUG;
		
		$self->{"config_data"}{$key} = $value;
		push @{ $self->{"file_fields"}{$file} }, $key;
		}
		
	close(CONFIG);
	
	$self->_validate_keys;
	
	return 1;
	}

=item parse_string( SCALAR_REF )

Parses the string inside the reference SCALAR_REF just as if
it found it in a file.

=cut

sub parse_string
	{
	my $self   = shift;
	my $string = shift;
	
	my @lines = split /\r?\n/, $$string;
	chomp( @lines );
	carp "A: Found " . @lines . " lines" if $DEBUG;
	
	while( my $line = shift @lines )
		{
		carp "1: Line is $line" if $DEBUG;

		CONT: {
		if ( $line =~ s/\\ \s* $//x )
			{
			carp "a: reading continuation line $lines[0]" if $DEBUG;
			$line .= shift @lines;
			carp "b: Line is $line" if $DEBUG;
			redo CONT unless @lines == 0;
			}
		}

		carp "2: Line is $line" if $DEBUG;
		
		chomp $line;
		next if $line =~ /^\s*(#|$)/; 

		carp "3: Line is $line" if $DEBUG;
				
		my ($key, $value) = &parse_line( $line );
		carp "Key:  '$key'   Value:  '$value'" if $DEBUG;
		
		$self->{"config_data"}{$key} = $value;
		}
			
	$self->_validate_keys;
	
	return 1;
	}
	
=item get( DIRECTIVE )

Returns the parsed value for that directive.  For directives
which did not have a value in the configuration file, C<get>
returns the empty string.

=cut

sub get 
	{
	my $self = shift;
	my $key  = shift;
	
	return $self->{"config_data"}{$key};
	}

=item set( DIRECTIVE, VALUE )

Sets the value for DIRECTIVE to VALUE.  The DIRECTIVE
need not already exist.  This overwrites previous 
values.

The VALUE must be a simple scalar.  It cannot be a reference.
If the VALUE is a reference, the function prints a warning
and returns false.

=cut

sub set 
	{
	my $self = shift;
	my( $key, $value ) = @_;
	
	if( ref $value )
		{
		$ERROR = "Second argument to set must be a simple scalar";
		if( $Warn )
			{
			carp $ERROR;
			return;
			}
		elsif( $Die )
			{
			croak $ERROR;
			}
		}
	
	$self->{"config_data"}{$key} = $value;
	}

=item unset( DIRECTIVE )

Remove the value from DIRECTIVE, which will still exist.  It's
value is undef.  If the DIRECTIVE does not exist, it will not
be created.  Returns FALSE if the DIRECTIVE does not already
exist, and TRUE otherwise.

=cut

sub unset
	{
	my $self = shift;
	my $key  = shift;
	
	return unless $self->exists( $key );
	
	$self->{"config_data"}{$key} = undef;
	
	return 1;
	}

=item remove( DIRECTIVE )

Remove the DIRECTIVE. Returns TRUE is DIRECTIVE existed
and FALSE otherwise.   

=cut

sub remove
	{
	my $self = shift;
	my $key  = shift;
	
	return unless $self->exists( $key );
	
	delete $self->{"config_data"}{$key};
	
	return 1;
	}

=item directives()

Returns a list of all of the directive names found in the configuration
file. The keys are sorted ASCII-betically.

=cut

sub directives
	{
	my $self = shift;

	my @keys = sort keys %{ $self->{"config_data"} };

	return @keys;
	}

=item exists( DIRECTIVE )

Return TRUE if the specified directive exists, and FALSE
otherwise.  

=cut

sub exists
	{
	my $self = shift;
	my $name = shift;
	
	return CORE::exists $self->{"config_data"}{ $name };
	}

=item clone

Return a copy of the object.  The new object is distinct
from the original so you can make changes to the new object
without affecting the old one.

=cut

# this is only the first stab at this -- from 35,000
# feet in coach class
# 
# I expect that the hash will be very simple.  Some keys
# might have a reference value, but that reference value
# will be "flat", so it won't have references in it.

sub clone
	{
	my $self = shift;
	
	my $clone = bless {}, ref $self;
	
	$clone->{"filenames"} = [ @{ $self->{"filenames"} } ];
	$clone->{"validkeys"} = [ @{ $self->{"validkeys"} } ];

	foreach my $file ( keys %{ $self->{"file_fields"} } )
		{
		$clone->{"file_fields"}{ $file } 
			= [ @{ $self->{"file_fields"}{ $file } } ];
		}
			
	foreach my $key ( $self->directives )
		{
		$clone->set( $key, $self->get( $key ) );
		}
				
	return $clone;
	}

=item save( FILENAME [ => ARRAY_REF [, FILENAME => ARRAY_REF ] ] );

The save method works in three ways, depending on the argument list.

With a single argument, the save function attempts to save all of the
field-value pairs of the object to the file named by the argument.

	$clone->save( "configrc" );

With two arguments, the method expects the second argument to be an
array reference which lists the directives to save in the file.

	$clone->save( "configrc" => [qw(Foo Bar)] );

With more than two arguments, the method expects filename-list pairs.
The method will save in each file the values in their respective 
array references.

	$clone->save(
		"configrc" => [qw(Foo Bar)],
		"global"   => [qw(Baz Quux)],
		);

In the last two cases, the method checks that the value for each pair
is an array reference before it affects any files.  It croaks if
any value is not an array reference.

Once the method starts writing files, it tries to write all of the
specified files. Even if it has a problem with one of them, it continues
onto the next one.  The method does not necessarily write the files
in the order they appear in the argument list, and it does not check
if you specified the same file twice.
	
=cut

sub save	
	{
	my $self = shift;
	my @args = @_;
	
	if( @args == 0 ) # no args!
		{
		carp "No arguments to method!";
		return;
		}
		
	if( @args == 1 )	# this is a single file
		{
		push @args, [ $self->directives ];
		}
		
	unless( @args % 2 == 0 ) { croak "Odd number of arguments" };
	
	my %hash = @args;
	
	foreach my $value ( values %hash )
		{
		croak "Argument is not an array reference"
			unless isa( $value, 'ARRAY' );
		}
		
	foreach my $file ( keys %hash )
		{
		carp $ERROR unless $self->_save( $file, $hash{$file} );
		}
		
	1;
	}
	
sub _save
	{
	my $self = shift;
	my $file = shift;
	
	my $directives = shift;

	unless( isa( $directives, 'ARRAY' ) )
		{
		$ERROR = 'Argument is not an array reference';
		return;
		}
		
	my $fh;
	
	unless( open $fh, "> $file" )
		{
		$ERROR = $!;
		
		return;
		}
		
	foreach my $directive ( @$directives )
		{
		print $fh ( join "\t", $directive, $self->get( $directive ) );
		print $fh "\n";
		}
		
	close $fh;
	
	return SUCCESS;
	}
	
# Internal methods

sub parse_line 
	{
	my $text = shift;
	
	my ($key, $value);
	
	# AWJ: Allow optional '=' or ' = ' between key and value:
	if ($text =~ /^\s*([^\s=]+)\s*[=]?\s*(['"]?)(.*?)\2\s*$/ ) 
		{
		( $key, $value ) = ( $1, $3 );
		} 
	else 
		{
		croak "Config: Can't parse line: $text\n";
		}
	
	return ($key, $value);
	}

sub _init_errors
	{
	%ERROR = ();
	$Error = undef;
	$ERROR = undef;
	}
	
# =item _validate_keys

# If any keys were declared when the object was constructed,
# check that those keys actually occur in the configuration file.
# This function croaks if a declared key does not exist.

# =cut

sub _validate_keys 
	{
	my $self = shift;
   
	if ( $self->{"validkeys"} )
		{
		my ($declared_key);
		my $declared_keys_ref = $self->{"validkeys"};

		foreach $declared_key ( @$declared_keys_ref )
			{
			unless ( $self->{"config_data"}{$declared_key} )
				{
				croak "Config: key '$declared_key' does not occur in file $self->{filename}\n";
      			}
         
         	carp "Key: $declared_key found.\n" if $DEBUG;
			}
		}

	return SUCCESS;
	}

=back

=head2 Package variables

=over 4

=item $Die

If set to a true value, all errors are fatal.

=item $ERROR

The last error message.

=item %ERROR

The error messages from unreadable files.  The key is
the filename and the value is the error message.

=item $Warn

If set to a true value, methods may output warnings.

=back

=head1 LIMITATIONS/BUGS

Directives are case-sensitive.

If a directive is repeated, the first instance will silently be
ignored.

=head1 CREDITS

Bek Oberin E<lt>gossamer@tertius.net.auE<gt> wote the original module

Kim Ryan E<lt>kimaryan@ozemail.com.auE<gt> adapted the module to make
declaring keys optional.  Thanks Kim.

Alan W. Jurgensen E<lt>jurgensen@berbee.comE<gt> added a change to allow
the NAME=VALUE format in the configuration file.

Andy Lester, E<lt>petdance@cpan.orgE<gt>, for maintaining the module
while brian was on active duty.

Adam Trickett, E<lt>atrickett@cpan.orgE<gt>, added multi-line support.
You might want to see his C<Config::Trivial> module.

Greg White has been a very patient user and tester.

=head1 SOURCE AVAILABILITY

This source is part of a SourceForge project which always has the
latest sources in CVS, as well as all of the previous releases.

	http://sourceforge.net/projects/brian-d-foy/
	
If, for some reason, I disappear from the world, one of the other
members of the project can shepherd this module appropriately.

=head1 AUTHORS

brian d foy, C<< <bdfoy@cpan.org> >>

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2002-2007 brian d foy.  All rights reserved.

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

1;
