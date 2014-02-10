# http.tcl
# Client-side HTTP for GET, POST, and HEAD commands.
# These routines can be used in untrusted code that uses the Safesock
# security policy.
# These procedures use a callback interface to avoid using vwait,
# which is not defined in the safe base.
#
# SCCS: @(#) http.tcl 1.2 97/01/22 13:12:02
#
# See the http.n man page for documentation

package provide http 1.0

if {[info commands "unsupported0"] == "unsupported0"} {
    rename unsupported0 copychannel
}
array set http {
    -proxyhost {}
    -proxyport {}
    -useragent {Tcl http client package 1.0}
    -proxyfilter httpProxyRequired
}
proc http_config {args} {
    global http
    if {[llength $args] == 0} {
	set result {}
	foreach name [lsort [array names http -*]] {
	    lappend result $name $http($name)
	}
	return $result
    } elseif {[llength $args] == 1} {
	set flag [lindex $args 0]
	if {[regexp -- {^-(proxyhost|proxyport|proxyfilter|agent)$} $flag]} {
	    return $http($flag)
	} else {
	    return -code error "Unknown option $flag, must be -proxyfilter, -proxyhost, -proxyport, or -useragent"
	}
    } else {
	foreach {flag value} $args {
	    switch -- $flag {
		-proxyhost -
		-proxyport -
		-proxyfilter -
		-useragent {
		    set http($flag) $value
		}
		default {
		    return -code error "Unknown option $flag, must be  -proxyfilter, -proxyhost, -proxyport, or -useragent"
		}
	    }
	}
    }
}

proc http_reset { token } {
    upvar #0 $token state
    set state(status) reset
    catch {fileevent $state(sock) readable {}}
    catch {eval $state(-command) {$token}}
    catch {close $state(sock)}
    catch {unset state}
}
proc http_get { url args } {
    global http
    if ![info exists http(uid)] {
	set http(uid) 0
    }
    set token http#[incr http(uid)]
    upvar #0 $token state
    http_reset $token
    array set state {
	-command 	{# }
	-blocksize 	8192
	-validate 	0
	-headers 	{}
	state		header
	meta		{}
	currentsize	0
	totalsize	0
        type            text/html
        body            {}
	status		""
    }
    foreach {flag value} $args {
	switch -- $flag {
	    -blocksize -
	    -channel -
	    -command -
	    -headers -
	    -progress -
	    -query -
	    -validate {
		set state($flag) $value
	    }
	    default {
		return -code error "Unknown option $flag: can be -blocksize, -channel, -command, -headers, -progress, -query, or -validate"
	    }
	}
    }
    if {! [regexp -nocase {^(http://)?([^/:]+)(:([0-9]+))?(/.*)} $url \
	    x proto host y port srvurl]} {
	error "Unsupported URL: $url"
    }
    if {[string length $port] == 0} {
	set port 80
    }
    if {[string length $proto] == 0} {
	set url http://$url
    }
    set state(url) $url
    if {![catch {$http(-proxyfilter) $host} proxy]} {
	set phost [lindex $proxy 0]
	set pport [lindex $proxy 1]
    }
    if {[info exists phost] && [string length $phost]} {
	set srvurl $url
	set s [socket $phost $pport]
    } else {
	set s [socket $host $port]
    }
    set state(sock) $s
    # Send data in cr-lf format, but accept any line terminators
    fconfigure $s -translation {auto crlf} -buffersize $state(-blocksize)
    # this is disallowed in safe interpreters, but the socket
    # is already in non-blocking mode in that case.
    catch {fconfigure $s -blocking off}
    set len 0
    set how GET
    if {[info exists state(-query)]} {
	set len [string length $state(-query)]
	if {$len > 0} {
	    set how POST
	}
    } elseif {$state(-validate)} {
	set how HEAD
    }
    puts $s "$how $srvurl HTTP/1.0"
    puts $s "Accept: */*"
    puts $s "Host: $host"
    puts $s "User-Agent: $http(-useragent)"
    foreach {key value} $state(-headers) {
	regsub -all \[\n\r\]  $value {} value
	set key [string trim $key]
	if {[string length $key]} {
	    puts $s "$key: $value"
	}
    }
    if {$len > 0} {
	puts $s "Content-Length: $len"
	puts $s "Content-Transfer-Encoding: x-url-encoding"
	puts $s ""
	fconfigure $s -translation {auto binary}
	puts $s $state(-query)
    } else {
	puts $s ""
    }
    flush $s
    fileevent $s readable [list httpEvent $token]
    return $token
}

 proc httpEvent {token} {
    upvar #0 $token state
    set s $state(sock)

    if [eof $s] then {
	close $s
	if {$state(state) == "header"} {
	    # Premature eof
	    set state(status) eof
	} else {
	    set state(status) ok
	}
	set state(state) eof
	eval $state(-command) {$token}
	return
    }
    if {$state(state) == "header"} {
	set n [gets $s line]
	if {$n == 0} {
	    set state(state) body
	    if ![regexp -nocase ^text $state(type)] {
		# Turn off conversions for non-text data
		fconfigure $s -translation binary
	    }
	} elseif {$n > 0} {
	    if [regexp -nocase {^content-type:(.+)$} $line x type] {
		set state(type) [string trim $type]
	    }
	    if [regexp -nocase {^content-length:(.+)$} $line x length] {
		set state(totalsize) [string trim $length]
	    }
	    if [regexp -nocase {^([^:]+):(.+)$} $line x key value] {
		lappend state(meta) $key $value
	    } elseif {[regexp ^HTTP $line]} {
		set state(http) $line
	    }
	}
    } else {
	if [catch {
#	    if [info exists state(-image)] {
#		$state(-image) config -channel $s
#	    } else
	    if {[info exists state(-channel)]} {
		set n [copychannel $s $state(-channel) $state(-blocksize)]
	    } else {
		set block [read $s $state(-blocksize)]
		set n [string length $block]
		if {$n >= 0} {
		    append state(body) $block
		}
	    }
	    if {$n >= 0} {
		incr state(currentsize) $n
	    }
	} err] {
	    set state(error) $err
	    http_reset $token
	} else {
	    if [info exists state(-progress)] {
		eval $state(-progress) {$token $state(totalsize) $state(currentsize)}
	    }
	}
    }
}
proc http_wait {token} {
    upvar #0 $token state
    if {![info exists state(status)] || [string length $state(status)] == 0} {
	vwait $token\(status)
    }
    return $state(status)
}

# Call http_formatQuery with an even number of arguments, where the first is
# a name, the second is a value, the third is another name, and so on.

proc http_formatQuery {args} {
    set result ""
    set sep ""
    foreach i $args {
	append result  $sep [httpMapReply $i]
	if {$sep != "="} {
	    set sep =
	} else {
	    set sep &
	}
    }
    return $result
}

# do x-www-urlencoded character mapping
# The spec says: "non-alphanumeric characters are replaced by '%HH'"
# 1 leave alphanumerics characters alone
# 2 Convert every other character to an array lookup
# 3 Escape constructs that are "special" to the tcl parser
# 4 "subst" the result, doing all the array substitutions
 
 proc httpMapReply {string} {
    global httpFormMap
    set alphanumeric	a-zA-Z0-9
    if ![info exists httpFormMap] {
	 
	for {set i 1} {$i <= 256} {incr i} {
	    set c [format %c $i]
	    if {![string match \[$alphanumeric\] $c]} {
		set httpFormMap($c) %[format %.2x $i]
	    }
	}
	# These are handled specially
	array set httpFormMap {
	    " " +   \n %0d%0a
	}
    }
    regsub -all \[^$alphanumeric\] $string {$httpFormMap(&)} string
    regsub -all \n $string {\\n} string
    regsub -all \t $string {\\t} string
    regsub -all {[][{})\\]\)} $string {\\&} string
    return [subst $string]
}

# Default proxy filter. 
 proc httpProxyRequired {host} {
    global http
    if {[info exists http(-proxyhost)] && [string length $http(-proxyhost)]} {
	if {![info exists http(-proxyport)] || ![string length $http(-proxyport)]} {
	    set http(-proxyport) 8080
	}
	return [list $http(-proxyhost) $http(-proxyport)]
    } else {
	return {}
    }
}
