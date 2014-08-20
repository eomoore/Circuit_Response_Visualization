$(document).ready(function(){

        $("#loading-div-background").hide();
        $("#loading-div-background").css({ opacity: 0.7 });

        $('#circ_form').submit(function(e){
                e.preventDefault();
                $("#content").hide();
                $("#loading-div-background").show();
                var form = $(this);
                var post_url = form.attr('action');
                var post_data = form.serialize();
                $.ajax({
                    type: 'POST',
                    url: post_url,
                    data: post_data,
                    success: function(msg) {
                        $('#content').html(msg);
                        $("#loading-div-background").hide();
                        $("#content").show();
                    }
                });
        });

        $('#annotate').click(function(){
                $("#content").hide();
                $("#loading-div-background").show();
                var form = $('#circ_form');
                var post_url = form.attr('action') + 'annotate/';
                // var post_url = form.attr('action');
                var post_data = form.serialize();
                $.ajax({
                    type: 'POST',
                    url: post_url,
                    data: post_data,
                    success: function(msg) {
                        $('#content').html(msg);
                        $("#loading-div-background").hide();
                        $("#content").show();
                    }
                });
        });
});
