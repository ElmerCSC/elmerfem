proc message { text { name message } } {

   toplevel .message
   place_window .message

   wm title .message $name

   frame .message.frame

   text .message.frame.text -width 42 -height 9 -font -Adobe-Helvetica-Medium-R-Normal-*-140-*

   button .message.frame.but -text "Close" -command "destroy .message"

   pack .message.frame
   pack .message.frame.text
   pack .message.frame.but -side bottom

   .message.frame.text insert end $text
}
