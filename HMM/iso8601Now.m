function timestampString = iso8601Now()

timestampString = strtrim(char(datetime('now','TimeZone','local','Format','uuuu-MM-dd'' ''HH:mm:ss XXX')));