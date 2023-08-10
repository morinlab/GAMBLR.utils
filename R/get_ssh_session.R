#' @title Get SSH Session.
#'
#' @description Create an ssh session to the GSC (requires active VPN connection)
#'
#' @details Using the ssh R package to create an ssh session.
#'
#' @param host Default is "gphost01.bcgsc.ca".
#'
#' @return An external pointer of class 'ssh_session'
#'
#' @import ssh
#' @export
#'
#' @examples
#' my_session = get_ssh_session()
#'
get_ssh_session = function(host="gphost01.bcgsc.ca"){

  if(!is.null(config::get("host"))){
    host = config::get("host")
  }

  if (!requireNamespace("ssh", quietly = TRUE)) {
    warning("The ssh package must be installed to use this functionality")
    #Either exit or do something that does not require ssh
    return(NULL)
  }
  message("you should also run this command to ensure the ssh library is loaded:\nlibrary(ssh)")
  session = ssh::ssh_connect(host=host)
  return(session)
}
