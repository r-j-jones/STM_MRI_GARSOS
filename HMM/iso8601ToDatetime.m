function t = iso8601ToDatetime(timestemapString)

t = datetime(timestemapString,'InputFormat','uuuu-MM-dd'' ''HH:mm:ss XXX','TimeZone','local');
t.Format = 'uuuu-MM-dd'' ''HH:mm:ss XXX';