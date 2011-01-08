//--------------------------------------------------
//
// Bedingte darstellung des SourceForge logos
//
function SFLogo()
{
  var x = getState();
  if (x==1)
  {
    document.write("<IMG src=\"sf.gif\" width=\"210\" height=\"62\" border=\"0\" alt=\"SourceForge.net Logo\" \/>");
  }
  else
  {
    document.write("<IMG src=\"http:\/\/sourceforge.net\/sflogo.php?group_id=137191&amp;type=5\" width=\"210\" height=\"62\" border=\"0\" alt=\"SourceForge.net Logo\" \/>");
  }
}

//--------------------------------------------------
//
// Statcounter code erzeugen
//
function statCounter()
{
  // statcounter variablen
  document.write("<script type=\"text\/javascript\" language=\"javascript\">");
  document.write("var sc_project=671117;");
  document.write("var sc_partition=5; ");
  document.write("var sc_security=\"f3fad66c\";");
  document.write("<\/script>");

  var sc_project=671117;
  var sc_partition=5;
  var sc_security="f3fad66c";

	
  // statcounter script	
  document.write("<script type=\"text\/javascript\" language=\"javascript\" src=\"http:\/\/www.statcounter.com\/counter\/counter.js\">");
  document.write("<\/script>");
}

//--------------------------------------------------
//
// abfragen ob der statistik blocking cookie aktiviert ist
//
function getState() 
{
  var val = "";
  if (document.cookie) 
  {
    var start = document.cookie.indexOf("=") + 1;
    var end   = document.cookie.indexOf(";");
    if (end == -1)
      end = document.cookie.length;

    val = document.cookie.substring(start, end);
  }

  return val;
}


//--------------------------------------------------
//
// Setzen des cookies
//
function setCookie (val) 
{
  var _expire = 86400 * 365 * 50000;   // 10 Jahre verfallszeit
  var _now = new Date();
  var _end = new Date(_now.getTime()+_expire);
  if (val==1)
  {	
    document.cookie = "nostat=1; expires=" + _end.toGMTString() + ";";
  }
  else
  {
    document.cookie = "nostat=0; expires=" + _end.toGMTString() + ";";
//    document.cookie = "nostat=0; expires=" + _end.toGMTString() + ";";
  }
}
