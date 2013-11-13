;+
; NAME:
;   mangadesign_version
; PURPOSE:
;   Return the version name for the mangadesign product
; CALLING SEQUENCE:
;   vers = mangadesign_version()
; OPTIONAL INPUTS:
;   /SIMPLE: Don't return svn revision number if version is trunk
; OUTPUTS:
;   vers       - Version name for the product mangadesign
; COMMENTS:
;   Depends on shell script in $MANGADESIGN_DIR/bin
;-
;------------------------------------------------------------------------------
function mangadesign_version,simple=simple
   cmd = "mangadesign_version"
   spawn, cmd, stdout, /noshell
   mangadesign_version = stdout[0]
   if keyword_set(simple) then mangadesign_version=(strsplit(mangadesign_version,/extract,' '))[0]
   return, mangadesign_version
end
;------------------------------------------------------------------------------
